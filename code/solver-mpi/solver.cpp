#include "solver.h"

/***************************************************************************************************
solverControl
****************************************************************************************************
Flow control to prepare and solve the system.
***************************************************************************************************/
void femSolver::solverControl(inputSettings* argSettings, triMesh* argMesh)
{
    int mype, npes;                // my processor rank and total number of processors
    MPI_Comm_rank(MPI_COMM_WORLD, &mype);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);
    if (mype==0) cout << endl << "=================== SOLUTION =====================" << endl;

    mesh = argMesh;
    settings = argSettings;
    int nec = mesh->getNec();

    for(int iec=0; iec<nec; iec++)
    {
        calculateJacobian(iec);
        calculateElementMatrices(iec);
    }
    applyDrichletBC();
    accumulateMass();
    explicitSolver();

    return;
}

/***************************************************************************************************
calculateJacobian
****************************************************************************************************
Compute and store the jacobian for each element.
***************************************************************************************************/
void femSolver::calculateJacobian(const int e)
{
    int myNode;       // node number for the current node
    double x[nen];    // x values for all the nodes of an element
    double y[nen];    // y values for all the nodes of an element

    triMasterElement* ME = mesh->getME(0);    // for easy access to the master element

    double J[2][2];              // Jacobian for the current element
    double detJ;                 // Jacobian determinant for the current element
    double invJ[2][2];           // inverse of Jacobian for the current element
    double dSdX[3];              // dSdx on a GQ point
    double dSdY[3];              // dSdy on a GQ point

    // collect element node coordinates in x[3] and y[3] matrices
    for (int i=0; i<nen; i++)
    {
        myNode =  mesh->getElem(e)->getLConn(i);
        x[i] = mesh->getLNode(myNode)->getX();
        y[i] = mesh->getLNode(myNode)->getY();
    }

    // for all GQ points detJ, dSDx[3] and dSdY[3] are determined.
    for (int p=0; p<nGQP; p++)
    {
        // Calculate Jacobian
        J[0][0] = ME[p].getDSdKsi(0)*x[0] + ME[p].getDSdKsi(1)*x[1] + ME[p].getDSdKsi(2)*x[2];
        J[0][1] = ME[p].getDSdKsi(0)*y[0] + ME[p].getDSdKsi(1)*y[1] + ME[p].getDSdKsi(2)*y[2];
        J[1][0] = ME[p].getDSdEta(0)*x[0] + ME[p].getDSdEta(1)*x[1] + ME[p].getDSdEta(2)*x[2];
        J[1][1] = ME[p].getDSdEta(0)*y[0] + ME[p].getDSdEta(1)*y[1] + ME[p].getDSdEta(2)*y[2];

        //Calculate determinant of Jacobian and store in mesh
        detJ = J[0][0]*J[1][1] - J[0][1]*J[1][0];
        mesh->getElem(e)->setDetJ(p, detJ);

        // Calculate inverse of Jacobian
        invJ[0][0] =  J[1][1]/detJ;
        invJ[0][1] = -J[0][1]/detJ;
        invJ[1][0] = -J[1][0]/detJ;
        invJ[1][1] =  J[0][0]/detJ;

        // Calculate dSdx and dSdy and store in mesh
        dSdX[0] = invJ[0][0]*ME[p].getDSdKsi(0) + invJ[0][1]*ME[p].getDSdEta(0);
        dSdX[1] = invJ[0][0]*ME[p].getDSdKsi(1) + invJ[0][1]*ME[p].getDSdEta(1);
        dSdX[2] = invJ[0][0]*ME[p].getDSdKsi(2) + invJ[0][1]*ME[p].getDSdEta(2);
        dSdY[0] = invJ[1][0]*ME[p].getDSdKsi(0) + invJ[1][1]*ME[p].getDSdEta(0);
        dSdY[1] = invJ[1][0]*ME[p].getDSdKsi(1) + invJ[1][1]*ME[p].getDSdEta(1);
        dSdY[2] = invJ[1][0]*ME[p].getDSdKsi(2) + invJ[1][1]*ME[p].getDSdEta(2);

        mesh->getElem(e)->setDSdX(p, 0, dSdX[0]);
        mesh->getElem(e)->setDSdX(p, 1, dSdX[1]);
        mesh->getElem(e)->setDSdX(p, 2, dSdX[2]);
        mesh->getElem(e)->setDSdY(p, 0, dSdY[0]);
        mesh->getElem(e)->setDSdY(p, 1, dSdY[1]);
        mesh->getElem(e)->setDSdY(p, 2, dSdY[2]);
    }

    return;
}

/***************************************************************************************************
void femSolver::calculateElementMatrices(const int e)
****************************************************************************************************
Compute the K, M and F matrices. Then accumulate the total mass into the node structure.
***************************************************************************************************/
void femSolver::calculateElementMatrices(const int e)
{
    int node;
    int D = settings->getD();
    double f = settings->getSource();
    double * xyz = mesh->getXyz();

    double totalM = 0.0;    // Total mass
    double totalDM = 0.0;   // Total diagonal mass
    double K[3][3];
    double M[3][3];
    double F[3];
    double radius, x, y;
    double rFlux = 0.01;
    double x_i[nen], y_i[nen];

    // First, fill M, K, F matrices with zero for the current element
    for(int i=0; i<nen; i++)
    {
        F[i] = 0.0;
        mesh->getElem(e)->setF(i, 0.0);
        mesh->getElem(e)->setM(i, 0.0);
        x_i[i] = xyz[mesh->getElem(e)-> getConn (i)*nsd + xsd ];
        y_i[i] = xyz[mesh->getElem(e)-> getConn (i)*nsd + ysd ];
        for(int j=0; j<nen; j++)
        {
            mesh->getElem(e)->setK(i, j, 0.0);
            K[i][j] = 0.0;
            M[i][j] = 0.0;
        }
    }

    // Now, calculate the M, K, F matrices
    for(int p=0; p<nGQP; p++)
    {
        const double zeta = mesh->getME(p)->getPoint(0); 
        const double eta  = mesh->getME(p)->getPoint(1);
        x = calculateRealPosition(zeta,eta,x_i);
        y = calculateRealPosition(zeta,eta,y_i);
        const double factor_F = mesh->getElem(e)->getDetJ(p) * mesh->getME(p)->getWeight() * calculateHeatSource(x,y);
        for(int i=0; i<nen; i++)
        {
            for(int j=0; j<nen; j++)
            {
                // Consistent mass matrix
                M[i][j] = M[i][j] +
                            mesh->getME(p)->getS(i) * mesh->getME(p)->getS(j) *
                            mesh->getElem(e)->getDetJ(p) * mesh->getME(p)->getWeight();
                // Stiffness matrix
                K[i][j] = K[i][j] +
                            D * mesh->getElem(e)->getDetJ(p) * mesh->getME(p)->getWeight() *
                            (mesh->getElem(e)->getDSdX(p,i) * mesh->getElem(e)->getDSdX(p,j) +
                            mesh->getElem(e)->getDSdY(p,i) * mesh->getElem(e)->getDSdY(p,j));
            }
        // Forcing matrix
        F[i] = F[i] + factor_F * mesh->getME(p)->getS(i);
        }
    }

    // For the explicit solution, it is necessary to have a diagonal mass matrix and for this,
    // lumping of the mass matrix is necessary. In order to lump the mass matrix, we first need to
    // calculate the total mass and the total diagonal mass:
    for(int i=0; i<nen; i++)
    {
        for(int j=0; j<nen; j++)
        {
            totalM = totalM + M[i][j];
            if (i==j)
                totalDM = totalDM + M[i][j];
        }
    }

    // Now the diagonal lumping can be done
    for(int i=0; i<nen; i++)
    {
        for(int j=0; j<nen; j++)
        {
            if (i==j)
                M[i][j] = M[i][j] * totalM / totalDM;
            else
                M[i][j] = 0.0;
        }
    }

    //Total mass at each node is accumulated on local node structure:
    for(int i=0; i<nen; i++)
    {
        node = mesh->getElem(e)->getLConn(i);
        mesh->getLNode(node)->addMass(M[i][i]);
    }

    // At this point we have the necessary K, M, F matrices as a member of femSolver object.
    // They must be hard copied to the corresponding triElement variables.
    for(int i=0; i<nen; i++)
    {
        node =  mesh->getElem(e)->getLConn(i);
        x = mesh->getLNode(node)->getX();
        y = mesh->getLNode(node)->getY();
        radius = sqrt(pow(x,2) + pow(y,2));
        if(radius < (rFlux - 1e-10))
        {
            mesh->getElem(e)->setF(i, F[i]);
        }
        mesh->getElem(e)->setM(i, M[i][i]);
        for(int j=0; j<nen; j++)
        {
            mesh->getElem(e)->setK(i,j,K[i][j]);
        }
    }

    return;
}

/***************************************************************************************************
void femSolver::applyDrichletBC()
****************************************************************************************************
if any of the boundary conditions set to Drichlet type
    visits all partition level nodes
        determines if any the nodes is on any of the side surfaces

***************************************************************************************************/
void femSolver::applyDrichletBC()
{
    int const nnc = mesh->getNnc();
    double * TG = mesh->getTG();
    double * xyz = mesh->getXyz();
    double x, y, radius;
    double alpha = settings->getAlpha();
    double beta = settings->getBeta();
    int mype;
    MPI_Comm_rank(MPI_COMM_WORLD, &mype);
    double rOuter = 0.1;
    double temp;
    this->nnSolved = 0;
    int partialnnSolved = 0;
    if (settings->getBC(1)->getType()==1)
    {
        for(int i=0; i<nnc; i++)
        {
            x = xyz[i*nsd+xsd];
            y = xyz[i*nsd+ysd];
            if (abs(x-0.0) <= 1E-10 || abs(x-1.0) <= 1E-10 ||abs(y-0.0) <= 1E-10 || abs(y-1.0) <= 1E-10)
            {
                    mesh->getNode(i)->setBCtype(1);
                    TG[i] = x*x + alpha * y*y + beta;
            }
            else
            {
                partialnnSolved += 1;
            }
        }
    }
    MPI_Allreduce(&partialnnSolved, &this->nnSolved, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    return;
}

/***************************************************************************************************
 * void femSolver::explicitSolver()
 ***************************************************************************************************
 *
 **************************************************************************************************/
void femSolver::explicitSolver()
{
    int const nn = mesh->getNn();
    int const nec = mesh->getNec();
    int const nnc = mesh->getNnc();
    int const nnl = mesh->getNnl();
    int mype, npes;                // my processor rank and total number of processors
    MPI_Comm_rank(MPI_COMM_WORLD, &mype);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);
    int const nIter = settings->getNIter();
    double const dT = settings->getDt();
    double T[3], MTnew[3];
    double mass;
    double * massG = mesh->getMassG();
    double * MTnewG = mesh->getMTnewG();
    double * MTnewL = mesh->getMTnewL();
    double * TG = mesh->getTG();
    double * TL = mesh->getTL();
    double massTmp, MTnewTmp;
    double MT;
    double Tnew;
    double partialL2error, globalL2error, initialL2error;
    double* M;                   // pointer to element mass matrix
    double* F;                   // pointer to element forcing vector
    double* K;                   // pointer to element stiffness matrix
    triElement* elem;            // temporary pointer to hold current element
    triNode* pNode;              // temporary pointer to hold partition nodes
    double maxT, minT, Tcur;
    minT = std::numeric_limits<double>::max();
    maxT = std::numeric_limits<double>::min();
    int nodeG[nen];
    double partialsumT = 0;
    double sumT;

    MPI_Win winMTnew;
    MPI_Win winTG;

    MPI_Win_create(MTnewG, nnc*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &winMTnew);
    MPI_Win_create(TG, nnc*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &winTG);
    MPI_Win_fence(0, winMTnew);
    MPI_Win_fence(0, winTG);

    localizeTemperature(winTG);

    for (int iter=0; iter<nIter; iter++)
    {

        // clear RHS MTnew at local node level
        for(int i=0; i<nnl; i++){
            MTnewL[i] = 0;
        }

        // Evaluate right hand side at element level
        for(int e=0; e<nec; e++)
        {
            elem = mesh->getElem(e);
            M = elem->getMptr();
            F = elem->getFptr();
            K = elem->getKptr();
            for(int i=0; i<nen; i++)
            {
                nodeG[i] = elem->getLConn(i);
                T[i] = TL[nodeG[i]];
            }

            MTnew[0] = M[0]*T[0] + dT*(F[0]-(K[0]*T[0]+K[1]*T[1]+K[2]*T[2]));
            MTnew[1] = M[1]*T[1] + dT*(F[1]-(K[3]*T[0]+K[4]*T[1]+K[5]*T[2]));
            MTnew[2] = M[2]*T[2] + dT*(F[2]-(K[6]*T[0]+K[7]*T[1]+K[8]*T[2]));

            // RHS is accumulated at local nodes
            MTnewL[nodeG[0]] += MTnew[0];
            MTnewL[nodeG[1]] += MTnew[1];
            MTnewL[nodeG[2]] += MTnew[2];
        }

        // local node level MTnew is transferred to partition node level (MTnewL to MTnewG)
        accumulateMTnew(winMTnew);

        // Evaluate the new temperature on each node on partition level
        partialL2error = 0.0;
        globalL2error = 0.0;
        for(int i=0; i<nnc; i++)
        {
            pNode = mesh->getNode(i);
            if(pNode->getBCtype() != 1)
            {
                massTmp = massG[i];
                MT = MTnewG[i];
                Tnew = MT/massTmp;
                partialL2error += pow(TG[i]-Tnew,2);
                TG[i] = Tnew;
                MTnewG[i] = 0;
            }
        }

        // Transfer new temperatures from partition to local node level (from TG to TL)
        localizeTemperature(winTG);

        // Reduce RMSerrors to node 0: does not work because only cpu 0 is going to break the loop
        MPI_Allreduce(&partialL2error, &globalL2error, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        globalL2error = sqrt(globalL2error/(double)this->nnSolved);

        if(iter==0)
        {
            initialL2error = globalL2error;
            if (mype==0){
                cout << "The initial error is: " << scientific << setprecision(15) << initialL2error << endl;
                cout << "Iter" << '\t' << "Time" << '\t' << "L2 Error" << '\t' << endl;
            }
        }
        globalL2error = globalL2error / initialL2error;

        if(mype==0 && iter%1000==0)
        {
            cout << "mype " << mype << " ";
            cout << iter << '\t';
            cout << fixed << setprecision(5) << iter*dT << '\t';
            cout << scientific << setprecision(15) << globalL2error << endl;
        }
    }

    MPI_Win_free(&winMTnew);
    MPI_Win_free(&winTG);

    partialsumT = 0.0;
    for (int i = 0; i < nnc; i++)
    {
        partialsumT += TG[i];
    }
    MPI_Allreduce(&partialsumT, &sumT, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    if (mype == 0) cout << "sumT " << sumT << endl;

    return;
}

/***************************************************************************************************
 * void femSolver::accumulateMass()
 ***************************************************************************************************
 *
 **************************************************************************************************/
void femSolver::accumulateMass()
{
    int n; //counter
    int in;     // node number in global
    int ipe;    // pe numnber where 'in' lies
    int inc;    // offset of in on ipe
    int npes;

    int const mnc = mesh->getMnc();
    int const nnl = mesh->getNnl();
    int const nnc = mesh->getNnc();
    int const nn = mesh->getNn();
    int const *  const nodeLToG = mesh->getNodeLToG();
    double * massG = mesh->getMassG();
    MPI_Comm_size(MPI_COMM_WORLD, &npes);

    double * buffer;
    buffer = new double [mnc];
    int nncTarget;

    MPI_Win win;

    MPI_Win_create(massG, nnc*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win);
    MPI_Win_fence(0, win);

    for(int ipes=0; ipes<npes; ipes++)
    {
        for (int i = 0; i < mnc; i++) buffer[i] = 0.0;
        for(int inl=0; inl<nnl; inl++)
            {
            in = nodeLToG[inl];
            ipe = in/mnc;
            inc = in%mnc;
            if (ipes == ipe) buffer[inc] = mesh->getLNode(inl)->getMass();;
            }
        if ((ipes+1)*mnc > nn) nncTarget = nn - ipes*mnc;
        else nncTarget = mnc;
        // cout << "ipes" << ipes << "npes" << npes << endl;
        // cout << "mype" << mype << " nncTarget " << nncTarget << endl;
        MPI_Accumulate(&buffer[0], nncTarget, MPI_DOUBLE, ipes, 0, nncTarget, MPI_DOUBLE, MPI_SUM, win);
        MPI_Win_fence(0, win);
    }

    MPI_Win_free(&win);

    delete[] buffer;

    return;
}

/***************************************************************************************************
 * void femSolver::accumulateMTnew()
 ***************************************************************************************************
 *
 **************************************************************************************************/
void femSolver::accumulateMTnew(MPI_Win win)
{
    int n; //counter
    int in;     // node number in global
    int ipe;    // pe numnber where 'in' lies
    int inc;    // offset of in on ipe
    int npes;

    int const mnc = mesh->getMnc();
    int const nnl = mesh->getNnl();
    int const nnc = mesh->getNnc();
    int const nn = mesh->getNn();
    int const *  const nodeLToG = mesh->getNodeLToG();
    MPI_Comm_size(MPI_COMM_WORLD, &npes);

    double * MTnewL = mesh->getMTnewL();
    double * MTnewG = mesh->getMTnewG();

    double * buffer;
    buffer = new double [mnc];
    int nncTarget;

    for(int ipes=0; ipes<npes; ipes++)
    {
        for (int i = 0; i < mnc; i++) buffer[i] = 0.0;

        for(int inl=0; inl<nnl; inl++)
        {
            in = nodeLToG[inl];
            ipe = in/mnc;
            inc = in%mnc;
            if (ipes == ipe) buffer[inc] = MTnewL[inl];
        }
        if ((ipes+1)*mnc > nn) nncTarget = nn - ipes*mnc;
        else nncTarget = mnc;
        MPI_Accumulate(&buffer[0], nncTarget, MPI_DOUBLE, ipes, 0, nncTarget, MPI_DOUBLE, MPI_SUM, win);
        MPI_Win_fence(0, win);

    }

    MPI_Win_fence(0, win);

    delete[] buffer;

    return;
}

/***************************************************************************************************
 * void femSolver::localizeTemperature()
 ***************************************************************************************************
 *
 **************************************************************************************************/
void femSolver::localizeTemperature(MPI_Win win)
{
    int n; //counter
    int in;     // node number in global
    int ipe;    // pe number where 'in' lies
    int inc;    // offset of in on ipe
    int npes;

    int const mnc = mesh->getMnc();
    int const nnl = mesh->getNnl();
    int const nnc = mesh->getNnc();
    int const nn = mesh->getNn();
    int const *  const nodeLToG = mesh->getNodeLToG();
    MPI_Comm_size(MPI_COMM_WORLD, &npes);


    double * TL = mesh->getTL();
    double * TG = mesh->getTG();

    double * buffer;
    buffer = new double [mnc];
    for (int i = 0; i < mnc; i++) buffer[i] = 0;
    int nncTarget;

    for(int ipes=0; ipes<npes; ipes++)
    {
        if ((ipes+1)*mnc > nn) nncTarget = nn - ipes*mnc;
        else nncTarget = mnc;
        MPI_Get(&buffer[0], nncTarget, MPI_DOUBLE, ipes, 0, nncTarget, MPI_DOUBLE, win);
        MPI_Win_fence(0, win);

        for(int inl=0; inl<nnl; inl++)
        {
            in = nodeLToG[inl];
            ipe = in/mnc;
            inc = in%mnc;
            if (ipes == ipe) TL[inl] = buffer[inc];
        }
    }

    MPI_Win_fence(0, win);

    delete[] buffer;

    return;
}
double femSolver:: calculateRealPosition(const double zeta, const double eta, const double* x) {
    return x[0] + (x[1] - x[0]) * zeta + (x[2] - x[0]) * eta;
}
double femSolver::calculateHeatSource(const double x, const double y) 
{
    double alpha = settings->getAlpha();
    double beta = settings->getBeta();

    return - 2 - 2*alpha;
}

/***************************************************************************************************
* void femSolver::setInitialTemperature()
***************************************************************************************************
*
**************************************************************************************************/
void femSolver::setInitialTemperature() 
{
    int const nnc = mesh->getNnc();
    double alpha = settings->getAlpha();
    double beta = settings->getBeta();
    double * TL = mesh->getTL();
    double * xyz = mesh->getXyz();
    triNode* pNode;
    double x, y;

    for(int i=0; i<nnc; i++)
    {
        y = xyz[i*nsd+ysd];
        pNode = mesh->getNode(i);
        if(pNode->getBCtype() != 1) 
        {
            TL[i] = 500;
        }
    }
    return;
}
