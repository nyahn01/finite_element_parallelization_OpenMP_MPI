#include "solver.h"
#include "settings.h"
/***************************************************************************************************
solverControl
****************************************************************************************************
Flow control to prepare and solve the system.
***************************************************************************************************/
void femSolver::solverControl(inputSettings* argSettings, triMesh* argMesh)
{
    cout << endl << "=================== SOLUTION =====================" << endl;

    mesh = argMesh;
    settings = argSettings;
    int ne = mesh->getNe();

    for(int iec=0; iec<ne; iec++)
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
    int myNode;     // node number for the current node
    double x[nen];  // x values for all the nodes of an element
    double y[nen];  // y values for all the nodes of an element

    triMasterElement* ME = mesh->getME(0);  // for easy access to the master element
    double * xyz = mesh->getXyz();

    double J[2][2];     // Jacobian for the current element
    double detJ;        // Jacobian determinant for the current element
    double invJ[2][2];  // inverse of Jacobian for the current element
    double dSdX[3];     // dSdx on a GQ point
    double dSdY[3];     // dSdy on a GQ point

    // collect element node coordinates in x[3] and y[3] matrices
    for (int i=0; i<nen; i++)
    {
        myNode =  mesh->getElem(e)->getConn(i);
        x[i] = xyz[myNode*nsd+xsd];
        y[i] = xyz[myNode*nsd+ysd];
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
    double alpha = settings->getAlpha();
    double beta  = settings->getBeta();
    double K[3][3];
    double M[3][3];
    double F[3];
    double radius;
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
        const double x = calculateRealPosition(zeta,eta,x_i);
        const double y = calculateRealPosition(zeta,eta,y_i);
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
        //F[i] = F[i] 
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
        node = mesh->getElem(e)->getConn(i);
        mesh->getNode(node)->addMass(M[i][i]);
    }

    // At this point we have the necessary K, M, F matrices as a member of femSolver object.
    // They must be hard copied to the corresponding triElement variables.
    for(int i=0; i<nen; i++)
    {
        node = mesh->getElem(e)->getConn(i);
        mesh->getElem(e)->setF(i, F[i]);
        mesh->getElem(e)->setM(i, M[i][i]);
        for(int j=0; j<nen; j++)
        {
            mesh->getElem(e)->setK(i,j,K[i][j]);
        }
    }

    if (e==3921)
    {
        for(int i=0; i<nen; i++)
        {
            for(int j=0; j<nen; j++)
            {
                cout << "K: " << K[i][j] << "\t";
            }
            cout << endl;
        }
        for(int i=0; i<nen; i++)
        {
            cout << "M: " << M[i][i] << endl;
        }
        for(int i=0; i<nen; i++)
        {
            cout << "F: " << F[i] << endl;
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
    int const nn = mesh->getNn();
    double * T = mesh->getT();
    double * xyz = mesh->getXyz();
    double x, y, radius;
    double alpha = settings->getAlpha();
    double beta  = settings->getBeta();
    double rOuter = 0.1;
    double temp;
    this->nnSolved = 0;
    //if any of the boundary conditions set to Drichlet BC
    if (settings->getBC(1)->getType()==1)
    {
        for(int i=0; i<nn; i++)
        {
            x = xyz[i*nsd+xsd];
            y = xyz[i*nsd+ysd];
            if (abs(x-0.0) <= 1E-10 || abs(x-1.0) <= 1E-10 ||abs(y-0.0) <= 1E-10 || abs(y-1.0) <= 1E-10)
            {
                    mesh->getNode(i)->setBCtype(1);
                    T[i] = x*x + alpha * y*y + beta;
            }
            else
            {
                this->nnSolved += 1;
            }
        }
    }
    cout << "nnSolved: " << this->nnSolved << endl;

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
    int const ne = mesh->getNe();
    int const nIter = settings->getNIter();
    double const dT = settings->getDt();
    double TL[3], MTnewL[3];
    double * massG = mesh->getMassG();
    double * MTnew = mesh->getMTnew();
    double * T = mesh->getT();
    double massTmp, MTnewTmp;
    double MT;
    double Tnew;
    double partialL2error, globalL2error, initialL2error;
    double* M;          // pointer to element mass matrix
    double* F;          // pointer to element forcing vector
    double* K;          // pointer to element stiffness matrix
    triElement* elem;   // temporary pointer to hold current element
    triNode* pNode;     // temporary pointer to hold partition nodes
    double maxT, minT, Tcur;
    minT = std::numeric_limits<double>::max();
    maxT = std::numeric_limits<double>::min();

    for (int iter=0; iter<nIter; iter++)
    {
        // clear RHS MTnew
        for(int i=0; i<nn; i++){
            MTnew[i] = 0;
        }
        if (iter==0)
        {

            setInitialTemperature();
        }

        // Evaluate right hand side at element level
        for(int e=0; e<ne; e++)
        {
            elem = mesh->getElem(e);
            M = elem->getMptr();
            F = elem->getFptr();
            K = elem->getKptr();
            for(int i=0; i<nen; i++)
            {
                TL[i] = T[elem->getConn(i)];
            }

            MTnewL[0] = M[0]*TL[0] + dT*(F[0]-(K[0]*TL[0]+K[1]*TL[1]+K[2]*TL[2]));
            MTnewL[1] = M[1]*TL[1] + dT*(F[1]-(K[3]*TL[0]+K[4]*TL[1]+K[5]*TL[2]));
            MTnewL[2] = M[2]*TL[2] + dT*(F[2]-(K[6]*TL[0]+K[7]*TL[1]+K[8]*TL[2]));

            // RHS is accumulated at local nodes
            MTnew[elem->getConn(0)] += MTnewL[0];
            MTnew[elem->getConn(1)] += MTnewL[1];
            MTnew[elem->getConn(2)] += MTnewL[2];
        }

        // Evaluate the new temperature on each node on partition level
        partialL2error = 0.0;
        globalL2error = 0.0;
        for(int i=0; i<nn; i++)
        {
            pNode = mesh->getNode(i);
            if(pNode->getBCtype() != 1)
            {
                massTmp = massG[i];
                MT = MTnew[i];
                Tnew = MT/massTmp;
                partialL2error += pow(T[i]-Tnew,2);
                T[i] = Tnew;
                MTnew[i] = 0;
            }
        }
        globalL2error = sqrt(partialL2error/this->nnSolved);

        if(iter==0)
        {
            initialL2error = globalL2error;
            cout << "The initial error is: " << initialL2error << endl;
            cout << "Iter" << '\t' << "Time" << '\t' << "L2 Error" << '\t' << endl;
        }

        globalL2error = globalL2error / initialL2error;

        if(iter%1000==0)
        {
            cout << iter << '\t';
            cout << fixed << setprecision(5) << iter*dT << '\t';
            cout << scientific << setprecision(5) << globalL2error << endl;
        }
    }
    return;
}


/***************************************************************************************************
* void femSolver::accumulateMass()
***************************************************************************************************
*
**************************************************************************************************/
void femSolver::accumulateMass()
{
    int nn = mesh->getNn();
    double * massG = mesh->getMassG();


    for(int i=0; i<nn; i++)
    {
        massG[i] = mesh->getNode(i)->getMass();
    }

    return;
}


/***************************************************************************************************
* void femSolver::calculateRealPosition()
***************************************************************************************************
*
**************************************************************************************************/
double femSolver::calculateRealPosition(const double zeta, const double eta, const double* x) 
{
    return x[0] + (x[1] - x[0]) * zeta + (x[2] - x[0]) * eta;
}

/***************************************************************************************************
* void femSolver::calculateHeatSource()
***************************************************************************************************
*
**************************************************************************************************/
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
    int const nn = mesh->getNn();
    double alpha = settings->getAlpha();
    double beta = settings->getBeta();
    double * T = mesh->getT();
    double * xyz = mesh->getXyz();
    triNode* pNode;
    double x, y, radius;

    for(int i=0; i<nn; i++)
    {
        y = xyz[i*nsd+ysd];
        pNode = mesh->getNode(i);
        if(pNode->getBCtype() != 1) 
        {
            //T[i] = x*x + alpha * y*y + beta;
            T[i] = 500;
        }
    }
    return;
}
