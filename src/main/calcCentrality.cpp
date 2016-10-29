#include "calcCentrality.hpp"


void getSubGraph(int nodes, std::vector<int> const * const network,
        std::vector<double> * const subGraph){
    arma::mat armaNet(nodes,nodes);
    armaNet.zeros();
    int i,j;
    for(i=0; i<nodes; i++)
        for(j=i; j<nodes; j++)
            if(network->at(i*nodes+j))
                armaNet(i,j)=armaNet(j,i)=1;

    arma::mat eigVecs;
    arma::colvec eigVals;

    arma::eig_sym(eigVals,eigVecs,armaNet);

    arma::mat eigVecsSq;
    arma::colvec eigValsExp;

    eigValsExp=arma::exp(eigVals);
    eigVecsSq=arma::pow(eigVecs,2.0);

    arma::colvec armaSubGraph;
    armaSubGraph=eigVecsSq*eigValsExp;

    subGraph->clear();
    for(i=0; i<nodes; i++)
        subGraph->push_back(armaSubGraph(i));
}


void getSubGraph(int nodes, std::vector<int> const * const network,
        std::vector<double> * const subGraph, int const deg){
    arma::sp_mat armaNet(nodes,nodes);
    int i,j;
    for(i=0; i<nodes; i++)
        for(j=i; j<nodes; j++)
            if(network->at(i*nodes + j))
                armaNet(i,j)=armaNet(j,i)=1;

    arma::sp_mat prod(nodes,nodes);
    prod.eye();

    subGraph->resize(nodes);
    std::fill(subGraph->begin(),subGraph->end(),0);
    for(i=0; i<deg; i++){
        prod*=armaNet/((double)(i+1));

        for(j=0; j<nodes; j++)
            subGraph->at(j)+=prod(j,j);
    }
}


// void getSubGraphEig(int nodes, std::vector<int> const * const network,
// 		    std::vector<double> * const subGraph, int const deg){
//     Eigen::SparseMatrix<double> armaNet(nodes,nodes);
//     int i,j;
//     for(i=0; i<nodes; i++)
//         for(j=i; j<nodes; j++)
//             if(network->at(i*nodes + j)){
//                 armaNet.insert(i,j)=1;
//                 if(i!=j)
//                     armaNet.insert(j,i)=1;
//             }

//     Eigen::SparseMatrix<double> prod(nodes,nodes);
//     Eigen::VectorXd diag;
//     prod.setIdentity();

//     subGraph->resize(nodes);
//     std::fill(subGraph->begin(),subGraph->end(),0);
//     for(i=0; i<deg; i++){
//         prod=prod*armaNet;
//         prod/=(double(i+1));

//         diag=prod.diagonal();
//         for(j=0; j<nodes; j++)
//             subGraph->at(j)+=diag(j);
//     }
// }




void getBetweenness(int nodes, std::vector<int> const * const network,
        std::vector<double> * const btwn){
    btwn->clear();
    btwn->resize(nodes);
    std::fill(btwn->begin(),btwn->end(),0.0);
    int s,v,w,lenPw,i;
    MyStack S;
    MyQueue Q;
    std::vector< std::vector<int> > P(nodes);
    std::vector<int> d(nodes);
    std::vector<double> delta(nodes),sigma(nodes);
    for(s=0; s<nodes; s++){

        S.clear();

        for(i=0; i<nodes; i++){
            P.at(i).clear();
            if(i!=s){
                sigma.at(i)=0;
                d.at(i)=-1;
            }
            else{
                sigma.at(i)=1;
                d.at(i)=0;
            }
        }

        Q.clear();
        Q.push(s);
        while(!Q.empty()){
            v=Q.pop();
            S.push(v);
            for(w=0; w<nodes; w++){
                if(network->at(w*nodes+v)){
                    // w found for the first time?
                    if(d.at(w) < 0){
                        Q.push(w);
                        d.at(w)=d.at(v)+1;
                    }
                    // shortest path to w via v?
                    if(d.at(w)==(d.at(v)+1)){
                        sigma.at(w)+=sigma.at(v);
                        P.at(w).push_back(v);
                    }
                }
            }
        }

        delta.clear();
        delta.resize(nodes);
        std::fill(delta.begin(),delta.end(),0.0);

        while(!S.empty()){
            w=S.pop();
            lenPw=P.at(w).size();
            for(i=0; i<lenPw; i++){
                v=P.at(w).at(i);
                delta.at(v)+=(sigma.at(v)/sigma.at(w))*(1.0+delta.at(w));
            }
            if(w!=s)
                btwn->at(w)+=delta.at(w);
        }
    }
}


int MyQueue::pop(){
    int val=vals.front();
    vals.pop_front();
    return val;
}

void MyQueue::push(int v){
    vals.push_back(v);
}

void MyQueue::clear(){
    vals.clear();
}

int MyQueue::empty(){
    return vals.empty();
}




int MyStack::pop(){
    int val=vals.front();
    vals.pop_front();
    return val;
}

void MyStack::push(int v){
    vals.push_front(v);
}

void MyStack::clear(){
    vals.clear();
}

int MyStack::empty(){
    return vals.empty();
}
