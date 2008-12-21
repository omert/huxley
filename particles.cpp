#include <iostream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <set>
#include "lapackpp.h"

using namespace std;

typedef LaGenMatDouble LaMatrix;
typedef LaSymmMatDouble LaSymMatrix;
typedef LaVectorDouble LaVector;

typedef vector<size_t> Permutation;

ostream& 
operator << (ostream& os, const Permutation& perm)
{
    for (size_t i = 0; i < perm.size(); ++i){
        if (i != 0)
            cout << " ";
        cout << perm[i];
    }
    return os;
}

//unordered pair
template<typename T>
struct Upair{
    Upair(const T& first, const T& second);
    bool operator < (const Upair<T>& other) const;
    T mFirst;
    T mSecond;
};

template<typename T>
bool 
Upair<T>::operator < (const Upair<T>& other) const
{
    if (mFirst == other.mFirst)
        return mSecond < other.mSecond;
    else
        return mFirst < other.mFirst;
}

template<typename T>
Upair<T>::Upair(const T& first, const T& second)
{
    if (first < second){
        mFirst = first;
        mSecond = second;
    }
    else{
        mFirst = second;
        mSecond = first;
    }
}


void
generatePermutations(const size_t n, const size_t k, 
                     Permutation perm,
                     vector<Permutation>& rPerms)
{
    if (perm.size() == k){
        rPerms.push_back(perm);
        return;
    }
    
    perm.push_back(0);
    for (size_t i = 0; i < n; ++i){
        bool used = false;
        for (size_t m = 0; m + 1 < perm.size(); ++m)
            if (perm[m] == i){
                used = true;
                break;
            }
        if (!used){
            perm.back() = i;
            generatePermutations(n, k, perm, rPerms);
        }
    }
}

void
generateCombinations(const size_t n, const size_t k, 
                     Permutation perm,
                     vector<Permutation>& rPerms)
{
    if (perm.size() == k){
        rPerms.push_back(perm);
        return;
    }
    
    perm.push_back(0);
    for (size_t i = 0; i < n; ++i){
        perm.back() = i;
        generateCombinations(n, k, perm, rPerms);
    }
}

set<Upair<size_t> >
permDiff(const Permutation& p1, const Permutation& p2)
{
    set<Upair<size_t> > ret;

    assert(p1.size() == p2.size());
    for (size_t i = 0; i < p1.size(); ++i)
        if (p1[i] != p2[i])
            ret.insert(Upair<size_t>(p1[i], p2[i]));

    return ret;
} 

void
printNice(ostream& os, const LaMatrix& M)
{
    for (int i = 0; i < M.rows(); ++i){
        for (int j = 0; j < M.cols(); ++j)
            os << setw(12) << setprecision(5) 
               << (fabs(M(i, j)) > 1e-10 ? M(i, j) : 0.0);
        os << endl;
    }
}

double
d(size_t i, size_t j)
{
    return i == j ? 1.0 : 0.0;
}

void
constructLargeMatrix(const LaSymMatrix& small, 
                     const vector<Permutation>& perms,
                     LaSymMatrix& rLarge)
{
    const double correlation = 1.0;
    size_t nPerm = perms.size();

    LaSymMatrix q = small;
    LaSymMatrix Q = rLarge;
    for (size_t ip1 = 0; ip1 < perms.size(); ++ip1)
        for (size_t ip2 = 0; ip2 < perms.size(); ++ip2){
            size_t i = perms[ip1][0];
            size_t j = perms[ip1][1];
            size_t k = perms[ip2][0];
            size_t l = perms[ip2][1];
            Q(ip1, ip2) = 
                q(j, l) * d(i, k) + 
                q(i, k) * d(j, l) +
                q(i, k) * d(i, l) * d(j, k) + 
                q(i, l) * d(i, k) * d(j, l);
        }
    cout << "Q:" << endl;
    printNice(cout, Q);
    

    for (size_t i = 0; i < perms.size(); ++i){
        for (size_t j = i + 1; j < perms.size(); ++j){
            set<Upair<size_t> > diff = permDiff(perms[i], perms[j]);
            if (diff.size() == 1)
                rLarge(i, j) = small(diff.begin()->mFirst, 
                                     diff.begin()->mSecond) * correlation;
            else
                rLarge(i, j) = 0.0;
//             if (perms[i][0] == perms[i][1] && perms[j][0] != perms[j][1])
//                 rLarge(i, j) = 0.0;
//             if (perms[j][0] == perms[j][1] && perms[i][0] != perms[i][1])
//                 rLarge(i, j) = 0.0;
//            cout << perms[i] << " -> " << perms[j] << ": " 
//                 << rLarge(i, j) << endl;
        }
    }
    for (size_t i = 0; i < perms.size(); ++i)
        for (size_t j = i; j < perms.size(); ++j){
            for (size_t k = 0; k < perms[i].size(); ++k)
                rLarge(i, j) += small(perms[i][k], perms[j][k]) * 
                    (1.0 - correlation);
        }
    for (size_t i = 0; i < nPerm; ++i){
        rLarge(i, i) = 0.0;
        for (size_t j = 0; j < nPerm; ++j)
            if (i != j)
                rLarge(i, i) -= rLarge(i, j);
    }
}

void
constructSmallMatrix(LaSymMatrix& rSmall)
{
    double scale = 1000.0;

    size_t nVert = rSmall.size(0);
    for (size_t i = 0; i < nVert; ++i)
        for (size_t j = 0; j < nVert; ++j)
            rSmall(i, j) = 0.0;

    for (size_t i = 0; i < nVert; ++i)
        for (size_t j = 0; j < i; ++j)
//            if (i - j == 1)
//                rSmall(i, j) = scale;
//            rSmall(i, j) = 1.0 * rand() / RAND_MAX > 0.5 ? scale : 0.0;
//            rSmall(i, j) = scale;
            rSmall(i, j) = scale * rand() / RAND_MAX;
            
            
    for (size_t i = 0; i < nVert; ++i)
        for (size_t j = 0; j < nVert; ++j)
            if (i != j)
                rSmall(i, i) -= rSmall(i, j);
} 

void
addSwitchParticleMatrix(const vector<Permutation>& perms, 
                size_t i1, size_t i2,
                double lambda,
                LaSymMatrix& rLarge)
{
    size_t nPerm = perms.size();
    LaMatrix switchMat = LaMatrix::zeros(nPerm, nPerm);
    for (size_t i = 0; i < nPerm; ++i){
        for (size_t j = i; j < nPerm; ++j)
            if (perms[i][i1] == perms[j][i2] && perms[i][i2] == perms[j][i1])
                switchMat(i, j) = switchMat(j, i) = 0.5;
        switchMat(i, i) -= 0.5;
    }

    LaMatrix comm = LaMatrix::zeros(nPerm, nPerm);
    comm = rLarge * switchMat - switchMat * rLarge;
    cout << "commutator:" << endl;
    printNice(cout, comm);

    for (size_t i = 0; i < nPerm; ++i)
        for (size_t j = i; j < nPerm; ++j)
            rLarge(i, j) += switchMat(i, j) * lambda;
}

class Switcher {
public:
    Switcher(size_t i1, size_t i2) : mi1(i1), mi2(i2) {}
    
    size_t operator()(size_t i){
        if (i == mi1)
            return mi2;
        else if (i == mi2)
            return mi1;
        else
            return i;
    }
    
private:
    const size_t mi1;
    const size_t mi2;
};

void
addSwitchVertexMatrix(const vector<Permutation>& perms, 
                      size_t i1, size_t i2,
                      double lambda,
                      LaSymMatrix& rLarge)
{
    Switcher switcher(i1, i2);

    size_t nPerm = perms.size();
    size_t nVert = perms[0].size();
    LaMatrix switchMat = LaMatrix::zeros(nPerm, nPerm);
    for (size_t i = 0; i < nPerm; ++i){
        Permutation switched(nVert);
        transform(perms[i].begin(), perms[i].end(), 
                  switched.begin(), switcher);
        for (size_t j = i; j < nPerm; ++j)
            if (switched == perms[j])
                switchMat(i, j) = switchMat(j, i) = 0.5;
        switchMat(i, i) -= 0.5;
    }

    LaMatrix comm = LaMatrix::zeros(nPerm, nPerm);
    comm = rLarge * switchMat - switchMat * rLarge;
    cout << "commutator:" << endl;
    printNice(cout, comm);

    for (size_t i = 0; i < nPerm; ++i)
        for (size_t j = i; j < nPerm; ++j)
            rLarge(i, j) += switchMat(i, j) * lambda;
}

void
showJointDistributions(const LaMatrix& largeEigenVectors, 
                       const vector<Permutation>& perms,
                       size_t nVert)
{
    cout << "Joint distributions" << endl;
    size_t nPerm = perms.size();
    for (size_t iV = 0; iV < nPerm; ++iV){
        LaMatrix dist = LaMatrix::zeros(nVert, nVert);
        for (size_t i = 0; i < nPerm; ++i)
            dist(perms[i][0], perms[i][1]) = largeEigenVectors(i, iV);
        cout << "Eigenvector " << iV << endl;
        printNice(cout, dist);
    }
}

void 
runtest(const size_t nVert, const size_t nPart)
{

    double rFact = 100000.0; //rounding factor

    LaSymMatrix small(nVert, nVert);
    constructSmallMatrix(small);

    vector<Permutation> perms;
    Permutation tempPerm;
    generatePermutations(nVert, nPart, tempPerm, perms);
//    generateCombinations(nVert, nPart, tempPerm, perms);
    size_t nPerm = perms.size();
    cout << nPerm << " permutations" << endl;

    cout << "permutations" << endl;
    for (size_t i = 0; i < perms.size(); ++i){
        cout << perms[i];
        cout << endl;
    }
    cout << endl;

    LaSymMatrix large(nPerm, nPerm);
    constructLargeMatrix(small, perms, large);
    
    cout << "Large Q:" << endl;
    printNice(cout, large);

//    addSwitchVertexMatrix(perms, 0, 1, 1.0, large);
//    addSwitchVertexMatrix(perms, 2, 3, 10.0, large);
//    addSwitchVertexMatrix(perms, 4, 5, 100.0, large);
    addSwitchParticleMatrix(perms, 0, 1, 10000.0, large);

    LaVector smallEigs(nVert);
    LaMatrix smallEigenVectors = small;
    LaEigSolve(small, smallEigs, smallEigenVectors);
    cout << "small Q:" << endl;
    printNice(cout, small);
    cout << endl << "small Q diagonalized:" << endl;
    printNice(cout, smallEigenVectors);
    cout << endl;
//    cout << smallEigs << endl;

    vector<double> smallEigsVec;
    for (size_t i = 0; i < nVert; ++i)
        smallEigsVec.push_back(round(-smallEigs(i) * rFact) / rFact);
    sort(smallEigsVec.begin(), smallEigsVec.end());
    
    cout << "small Q eigenvalues:" << endl;
    for (size_t i = 0; i < nVert; ++i)
        cout << smallEigsVec[i] << " ";
    cout << endl << endl << "sums of Q eigenvalues:" << endl;
    for (size_t i = 0; i < nVert; ++i){
        for (size_t j = 0; j <= i; ++j)
            cout << setw(12) << setprecision(6) 
                 << smallEigsVec[i] + smallEigsVec[j] << " ";
        cout << endl;
    }
    cout << endl;

    LaVector largeEigs(nPerm);
    LaMatrix largeEigenVectors = large;
    LaEigSolve(large, largeEigs, largeEigenVectors);
    LaMatrix largeEigsTransposed(1, nPerm);
    for (size_t i = 0; i < nPerm; ++i)
        largeEigsTransposed(0, i) = round(largeEigs(i) * rFact) / rFact;
    printNice(cout, largeEigsTransposed);
    cout << endl;
    printNice(cout, largeEigenVectors);
//    cout << largeEigs << endl;

    vector<double> largeEigsVec;
    for (size_t i = 0; i < nPerm; ++i)
        largeEigsVec.push_back(round(largeEigs(i) * rFact) / rFact);
    sort(largeEigsVec.begin(), largeEigsVec.end());

    cout << "big Q eigenvalues:" << endl;
    double lastEig = largeEigsVec[0];
    double eigCount = 1;
    for (size_t i = 1; i < nPerm; ++i)
        if (largeEigsVec[i] != lastEig){
            cout << lastEig << "(" << eigCount << ") ";
            eigCount = 1;
            lastEig = largeEigsVec[i];
        }
        else
            ++eigCount;
    cout << lastEig << "(" << eigCount << ") ";
    cout << endl << endl;

//    cout << "big Q eigenvectors:" << endl;    
//    printNice(cout, largeEigenVectors);
//   cout << endl;

    cout << "marginal distributions:" << endl;
    for (size_t iPart = 0; iPart < nPart; ++iPart){
        LaMatrix marginal = LaMatrix::zeros(nVert, nPerm);
        for (size_t iEig = 0; iEig < nPerm; ++iEig){
            for (size_t iPerm = 0; iPerm < nPerm; ++iPerm)
                marginal(perms[iPerm][iPart], iEig) += 
                    largeEigenVectors(iPerm, iEig);
            double invNorm = 0.0;
            for (size_t iVert = 0; iVert < nVert; ++iVert)
                invNorm += pow(marginal(iVert, iEig), 2);
            if (invNorm > 1e-4)
                invNorm = 1.0 / sqrt(invNorm);
            else
                invNorm = 0.0;
            for (size_t iVert = 0; iVert < nVert; ++iVert)
                marginal(iVert, iEig) *= invNorm;
        }
        cout << "particle " << iPart << endl;
        printNice(cout, marginal);
        cout << endl;
    }

    showJointDistributions(largeEigenVectors, perms, nVert);

    LaMatrix testVecs(nPerm, nPerm);
    for (size_t iEigComb = 0; iEigComb < nPerm; ++iEigComb){
        LaVector f(nVert, 1);
        for (size_t iv = 0; iv < nVert; ++iv)
            f(iv) = smallEigenVectors(iv, perms[iEigComb][0]);
        LaVector g(nVert, 1);
        for (size_t iv = 0; iv < nVert; ++iv)
            g(iv) = smallEigenVectors(iv, perms[iEigComb][1]);
        for (size_t iPerm = 0; iPerm < nPerm; ++iPerm){
            size_t i = perms[iPerm][0];
            size_t j = perms[iPerm][1];
            if (perms[iEigComb][1] > perms[iEigComb][0])  //symmetric
                testVecs(iEigComb, iPerm) = 
                    f(i) * g(j) + f(j) * g(i) + 
                    2.0 / (nVert - 2.0) * 
                    (f(i) * g(i) + f(j) * g(j));
            else  //anti-symmetric
                testVecs(iEigComb, iPerm) = f(i) * g(j) - f(j) * g(i);
        }
    }

    for (size_t i = 0;  i < nPerm; ++i){
        double invNorm = 0.0;
        for (size_t j = 0;  j < nPerm; ++j)
            invNorm += testVecs(i, j) * testVecs(i, j);
        if (invNorm > 1e-4)
            invNorm = 1.0 / sqrt(invNorm);
        else
            invNorm = 0.0;
        for (size_t j = 0;  j < nPerm; ++j)
            testVecs(i, j) *= invNorm;
    }

    cout << endl << "test vecs" << endl;;
    printNice(cout, testVecs);


    testVecs = testVecs * largeEigenVectors;

    for (size_t i = 0;  i < nPerm; ++i)
        for (size_t j = 0;  j < nPerm; ++j)
            if (fabs(testVecs(i, j)) < 1e-4)
                testVecs(i, j) = 0.0;

    cout << endl << "test vecs" << endl;;
    printNice(cout, testVecs);
    

    double e1 = smallEigsVec[1];
    double e2 = largeEigsVec[1];
    if (fabs(e1 - e2) != 0.0){
        cout << endl << largeEigsVec[1] << " " << smallEigsVec[1] << endl
             << small << endl;
        cout << endl; 
        sleep(2);
        fflush(NULL);
//        exit(0);
    }
    cout << "----------------------------------------------" << endl;
    fflush(NULL);
}

int main(int argc, char* argv[])
{
    size_t nVert = 4; //number of vertices in graph
    size_t nPart = 2; //number of particles
    
    if (argc > 1)
        nVert = atoi(argv[1]);
    if (argc > 2)
        nPart = atoi(argv[2]);

    cout << endl << "=========================================" << endl
         << nVert << " vertices, " << nPart << " particles" << endl;
    if (nPart < 1 || nVert < nPart){
        cout << "must have at least one particle, and at least as many "
             << "vertices as particles" << endl;
        fflush(0);
        exit(-1);
    }

    for(size_t i = 0; i < 100; ++i)
        runtest(nVert, nPart);
    return 0;
}
