#include <iostream>
#include <fstream>
#include <string.h>
using namespace std;
int viz[101];
class NFA
{
    int q;
    int f;
    int *F;
    int v;
    char *V;
    int g;
    struct tranz
    {
        int a,b;
        char c;
    } *G;
    int sc;
    int exista;
public:
    NFA(char *fisier)
    {
        int i;
        ifstream in(fisier);
        in>>q>>f;
        F=new int[f];
        for (i=0; i<f; i++)
            in>>F[i];
        in>>v;
        V=new char[v];
        for (i=0; i<v; i++)
            in>>V[i];
        in>>g;
        G=new tranz[g];
        for (i=0; i<g; i++)
        {
            in>>G[i].a>>G[i].c>>G[i].b;
        }
        sc=0;
        exista=0;
    }

    int GetG()
    {
        return g;
    }
    int GetQ()
    {
        return q;
    }
    int GetGa(int i)
    {
        return G[i].a;
    }
    int GetGb(int i)
    {
        return G[i].b;
    }
    char GetGc(int i)
    {
        return G[i].c;
    }
    char GetV(int i)
    {
        return V[i];
    }
    int GetV()
    {
        return v;
    }
    int GetF()
    {
        return f;
    }
    int GetFvector(int i)
    {
        return F[i];
    }
    void SetF(int x)
    {
        f = x;
    }
    void SetFvector(int x[101], int i)
    {
        F[i] = x[i];
    }
    void SetQ(int x)
    {
        q = x;
    }
    void SetG(int x)
    {
        g = x;
    }
    void SetTriplete(int x, char z, int y, int i)
    {
        G[i].a = x;
        G[i].b = y;
        G[i].c = z;
    }
    void SetVvector(char c, int i)
    {
        V[i] = c;
    }
    void SetV(int i)
    {
        v = i;
    }
    void show();
    bool verificare(char *cuvant, int sc, int pozitia_in_cuvant);
};

void NFA::show()
{
    int i;
    for (i=0; i<g; i++)
        cout<<G[i].a<<" "<<G[i].c<<" "<<G[i].b<<"\n";
}
bool NFA::verificare(char *cuvant, int sc, int pozitia_in_cuvant)
{
    int i;
    if ((pozitia_in_cuvant == strlen(cuvant)) && (sc <= q))
    {
        for (i=0; i<f; i++)
            if (sc==F[i])
            {
                exista=1;
                return true;
            }
    }
    else
    {
        for (i=0; i<g; i++)
        {
            if (exista==1)
                return true;
            if ((G[i].a==sc)&&(G[i].c==cuvant[pozitia_in_cuvant]))
            {
                verificare(cuvant,G[i].b,pozitia_in_cuvant+1);
            }
            else if(G[i].a==sc && G[i].c=='$')
                verificare(cuvant, G[i].b, pozitia_in_cuvant);
        }
        return false;
    }
}
void DF$(int k, NFA x, int n)
{
    int i;
    viz[k] = 1;
    for(i=0; i<n; i++)
        if(viz[i] == 0)
            for(int j=0; j<x.GetG(); j++)
                if(x.GetGa(j) == k && x.GetGb(j) == i && x.GetGc(j) == '$')
                    DF$(i, x, n);
}
void DFcaracter(int k, NFA x, int n, char caracter)
{
    int i, j;
    viz[k]++;
    for(i=0; i<n; i++)
        if(viz[i] <= 0)
            for(j=0; j<x.GetG(); j++)
                if(x.GetGa(j) == k && x.GetGb(j) == i)
                {
                    if(x.GetGc(j) == '$')
                    {
                        viz[i]--;
                        DFcaracter(i, x, n, caracter);
                    }
                    else if(x.GetGc(j) == caracter)
                    {
                        viz[i]++;
                    }
                }
}

NFA lnfa_to_nfa(NFA automat){
    int n = automat.GetQ(), fr1[n][n], fr2[n][n], costerg, vecaux[n], fr[n][n], p, i1, inlocuitor[n], reun[n], m1, m2, ok1, ok2, c, lin, caux, cstaux, csterg, i, j, l, m, staux[n], ok, okj, okl, f = automat.GetF(), F[101], VecMat[automat.GetV()][n][n], aux, sterg[n];
    for(i=0; i<f; i++)
        F[i] = automat.GetFvector(i);
    for(i=0; i<n; i++)          ///initializarea matricei lambda-inchidere
        for(j=0; j<n; j++)
            fr1[i][j] = -1;
    for(i=0; i<n; i++)          ///calcularea lamda-inchiderii
    {
        for(j=0; j<n; j++)
            viz[j] = 0;
        DF$(i, automat, n);
        c = 0;
        for(j=0; j<n; j++)
            if(viz[j] > 0)
            {
                fr1[i][c] = j;
                c++;
            }
    }

    for(i=0; i<n; i++)                      ///initializarea vectorului ce va determia noile stari finale
        reun[i] = 0;
    for(i=0; i<n; i++)
        for(j=0; j<n; j++)
            if(fr1[i][j] != -1)
                for(l=0; l<f; l++)
                    if(fr1[i][j] == F[l])
                        reun[i] = 1;
    c = 0;
    for(i=0; i<n; i++)                      ///memorarea noilor stari fnale
        if(reun[i] == 1)
        {
            F[c] = i;
            c++;
        }
    f = c;
    for(l=0; l<automat.GetV(); l++)
        if(automat.GetV(l) != '$')
        {
            for(i=0; i<n; i++)              ///initializarea caracter-inchiderii
                for(j=0; j<n; j++)
                    fr2[i][j] = -1;
            for(i=0; i<n; i++)
            {
                for(j=0; j<n; j++)
                    if(j != i)
                        viz[j] = 0;
                    else viz[i] = -1;
                DFcaracter(i, automat, n, automat.GetV(l));
                c = 0;
                for(j=0; j<n; j++)
                    if(viz[j] > 0)
                    {
                        fr2[i][c] = j;
                        c++;
                    }
            }
            for(i=0; i<n; i++)
            {
                reun[i] = 0;
            }
            for(i=0; i<n; i++)
                for(j=0; j<n; j++)
                    fr[i][j] = -1;
            for(j=0; j<n; j++)
            {
                i = 0;
                while(fr2[j][i] != -1)
                {
                    m = 0;
                    while(fr1[fr2[j][i]][m] != -1)
                    {
                        reun[fr1[fr2[j][i]][m]] = 1;
                        m++;
                    }
                    i++;
                }
                c = 0;
                for(m=0; m<n; m++)
                {
                    if(reun[m] != 0)
                    {
                        fr[j][c] = m;
                        c++;
                    }
                }
                for(m=0; m<n; m++)
                {
                    reun[m] = 0;
                }
            }
            for(i=0; i<n; i++)
            {
                for(j=0; j<n; j++)
                {
                    VecMat[l][i][j] = fr[i][j];
                }

            }
        }
    for (i=0; i<n; i++)
        inlocuitor[i] = -1;
    for (i=0; i<n; i++)
        sterg[i] = i;
    costerg = n;
    for(i=0; i<automat.GetV(); i++)
    {
        if(automat.GetV(i) != '$'){
            c = 0;
            for(j=0; j<n-1; j++)
            {
                for(l=j+1; l<n; l++)
                {
                    ok = 1;
                    for(m=0; m<n; m++)
                        if(VecMat[i][j][m] != VecMat[i][l][m])
                            ok = 0;
                    if(ok == 1)
                    {
                        okj = 0;
                        okl = 0;

                        for(m=0; m<f; m++)
                        {
                            if(F[m] == j)
                                okj = 1;
                            if(F[m] == l)
                                okl = 1;
                        }
                        if(okj == okl)
                        {
                            for(i1=0; i1<i; i1++)
                                for(m=0; m<n; m++)
                                    if(VecMat[i1][j][m] != VecMat[i1][l][m])
                                        ok = 0;
                            if(ok == 1)
                            {
                                vecaux[c] = l;
                                inlocuitor[l] = j;
                                c++;
                            }

                        }

                    }

                }

            }
        }

        for(j=0; j<c-1; j++)
            for(l=j+1; l<c; l++)
                if(vecaux[j]>vecaux[l])
                {
                    aux = vecaux[j];
                    vecaux[j] = vecaux[l];
                    vecaux[l] = aux;
                }
        caux = 0;
        csterg = 0;
        cstaux = 0;
        while(caux < c && csterg < costerg)
        {

            if(sterg[csterg] == vecaux[caux])
            {
                staux[cstaux] = sterg[csterg];
                caux++;
                csterg++;
                cstaux++;
            }
            else if(sterg[csterg] < vecaux[caux])
                csterg++;
            else caux++;
        }
        for(m=0; m<cstaux; m++)
            sterg[m] = staux[m];
        costerg = cstaux;
        for(m=0; m<n; m++)
        {
            if(inlocuitor[m] != -1)
            {
                aux = inlocuitor[m];
                inlocuitor[m] = -1;
                for(l=0; l<cstaux; l++)
                {
                    if(m == sterg[l])
                        inlocuitor[m] = aux;
                }
            }

        }
    }
    for(lin=costerg-1; lin>=0; lin--)
    {
        //cout<<"Vom sterge linia cu nr. "<<sterg[lin]<<endl;
        for(i=0; i<automat.GetV(); i++)
            if(automat.GetV(i) != '$')
            {
                for(l=0; l<n; l++)
                {
                    VecMat[i][sterg[lin]][l] = -2;
                }
                for(j=0; j<n; j++)
                {
                    ok1 = -1;
                    ok2 = 0;
                    for(m=0; m<n; m++)
                        {
                            if(VecMat[i][j][m] == sterg[lin])
                            {
                                ok1 = m;
                            }
                            if(VecMat[i][j][m] == inlocuitor[sterg[lin]])
                            {
                                ok2 = 1;
                            }
                        }
                        if(ok1 != -1)
                        {
                            if(ok2 == 1)
                            {
                                for(l=ok1; l<n-1; l++)
                                    VecMat[i][j][l] = VecMat[i][j][l+1];
                                VecMat[i][j][n-1] = -1;
                            }
                            else
                            {
                                VecMat[i][j][ok1] = inlocuitor[sterg[lin]];
                                for(m1=0; m1<n-1; m1++)
                                    for(m2=m1+1; m2<n; m2++)
                                        if(VecMat[i][j][m2] != -1)
                                            if(VecMat[i][j][m1] > VecMat[i][j][m2])
                                            {
                                                aux = VecMat[i][j][m1];
                                                VecMat[i][j][m1] = VecMat[i][j][m2];
                                                VecMat[i][j][m2] = aux;
                                            }
                            }
                        }

                    }
                /* ///afisarea matricilor cu liniile corespunzatoare sterse
                for(int j=0; j<n; j++)
                {
                    if(VecMat[i][j][0] != -2)
                    {
                        cout<<j<<": ";
                        for(int l=0; l<n; l++)
                            cout<<VecMat[i][j][l]<<" ";
                        cout<<endl;
                    }
                }
                cout<<endl;
                */
            }
    }
    automat.SetQ(n);

    p = 0;
    while(p < f){
        if(VecMat[0][F[p]][0] == -2){
            for(i=p; i<f-1; i++)
                F[i] = F[i + 1];
            f--;
        }else p++;
    }

    automat.SetF(f);
    for(i=0; i<f; i++)
        automat.SetFvector(F, i);
    p = 0;
    for(l=0; l<automat.GetV(); l++)
        if(automat.GetV(l) != '$')
        {
            for(i=0; i<n; i++)
                for(j=0; j<n; j++)
                    if(VecMat[l][i][j] > -1)
                    {
                        automat.SetTriplete(i, automat.GetV(l), VecMat[l][i][j], p);
                        p++;
                    }
        }
    automat.SetG(p);

    for(i=0; i<automat.GetV()-1; i++){
        if(automat.GetV(i) == '$'){
            for(j=i; j<automat.GetV(); j++){
                automat.SetVvector(automat.GetV(i+1), i);
            }
        }
    }
    automat.SetVvector(NULL, automat.GetV()-1);
    automat.SetV(automat.GetV()-1);
    /* ///afisare de ajutor
    for(i=0; i<automat.GetV()-1; i++)
    {
        cout<<automat.GetV(i)<<": "<<endl;
        for(int j=0; j<n; j++)
        {
            if(VecMat[i][j][0] != -2)
            {
                cout<<j<<": ";
                for(int l=0; l<n; l++)
                    cout<<VecMat[i][j][l]<<" ";
                cout<<endl;
            }
        }
        cout<<endl;
    }
    */
    return automat;
}


NFA nfa_to_dfa(NFA automat){
    int coada[101][101], j, k, l, v[100], nr, leg[101][3], i, u, p, c, n = automat.GetQ(), ok, VecMat[50][50][50], f, F[100];   ///cel mai nefavorabil caz este ca acest sir sa aiba 2^n elemente
    f = automat.GetF();
    for(i=0; i<f; i++){
        F[i] = automat.GetFvector(i);
    }
    for(i=0; i<automat.GetV(); i++)
    {
        for(int j=0; j<n; j++)
        {
            for(int l=0; l<n; l++)
                VecMat[i][j][l] = -1;
        }
    }
    for(l=0; l<automat.GetG(); l++){
        for(i=0; i<automat.GetV(); i++)
            {
            if(automat.GetV(i) == automat.GetGc(l)){
                p = 0;
                while(VecMat[i][automat.GetGa(l)][p] != -1){
                    p++;
                }
                VecMat[i][automat.GetGa(l)][p] = automat.GetGb(l);
            }
        }
    }

    coada[0][0] = 1;     ///primul element al cozii este mereu 0 deoarece acesta nu va fi niciodata sters
    c = 0;
    p = 0;
    u = 0;
    do{
        for(k=0; k<automat.GetV(); k++){
            ///prelucrarea lui coada[p]
            for (i=0; i<n; i++)
                v[i]=0;
            ok = 0;
            for (i=0; i<n; i++)
                if (coada[p][i]==1)
                {
                    for (j=0; j<n; j++)
                        if (VecMat[k][i][j]!= -1)
                        {
                            v[VecMat[k][i][j]] = 1;
                            ok = 1;
                        }
                }
            if (ok)
            {

                for (i=0; i<=u; i++)
                {
                    nr=0;
                    for (j=0; j<n; j++)
                        if (coada[i][j]==v[j]) nr++;
                    if (nr==n)
                    {
                        ok=0;
                        leg[c][0] = p;
                        leg[c][1] = k;
                        leg[c][2] = i;
                        c++;
                    }
                }
                if (ok)
                {
                    u++;
                    leg[c][0] = p;
                    leg[c][1] = k;
                    leg[c][2] = u;
                    c++;
                    for (j=0; j<n; j++)
                        coada[u][j]=v[j];
                }
            }
        }
        p++;
    }while(p <= u);
    n = u + 1;

    p = 0;
    for(i=0; i<n; i++){
        for(j=0; j<n; j++){
            if(coada[i][j] == 1){
                ok = 0;
                for(l=0; l<f; l++)
                    if(F[l] == j){
                        ok = 1;
                    }

                if(ok == 1){
                    for(l=0; l<p; l++)
                        if(v[l] == i){
                            ok = 0;
                        }

                }
                if(ok == 1){
                    v[p] = i;
                    p++;
                }
            }
        }
    }
    f = p;
    for(i=0; i<f; i++)
        F[i] = v[i];

    for(i=0; i<c; i++)
    {
        v[i] = i;
    }
    /* ///afisare de ajutor
    p = 0;
    while(p <= u)
    {
        for (j=0; j<n; j++)
            if (coada[p][j]!=0)
                cout<<j;
        cout<<" ";
        p++;
    }
    */
    automat.SetQ(n);
    automat.SetF(f);
    for(i=0; i<f; i++)
        automat.SetFvector(F, i);
    automat.SetG(c);
    for(i=0; i<c; i++)
    {
        automat.SetTriplete(leg[i][0], automat.GetV(leg[i][1]), leg[i][2], i);
    }
    return automat;
}


NFA dfa_to_dfamin(NFA automat){
    int n = automat.GetQ(), nr, si, Mat[n][automat.GetV()], i, j, l, fr[n], multimi[n+1][n+1], VecMat[automat.GetV()][n][n], p, u, k, m, frecv[n+1];
    bool MatFin[n][n], ok1, ok2;
    for(i=0; i<n; i++)
        for(j=0; j<automat.GetV(); j++)
            if(automat.GetV(j) != '$')
                Mat[i][j] = -1;
    for(i=0; i<automat.GetG(); i++){
        for(j=0; j<automat.GetV(); j++)
            if(automat.GetV(j) == automat.GetGc(i)){
                Mat[automat.GetGa(i)][j] = automat.GetGb(i);
            }
    }
    for(i=0; i<n; i++){
        for(j=0; j<i; j++){
            ok1 = false;
            ok2 = false;
            for(l=0; l<automat.GetF(); l++){
                if(i == automat.GetFvector(l))
                    ok1 = true;
                if(j == automat.GetFvector(l))
                    ok2 = true;
            }
            if(ok1 == ok2){
                MatFin[i][j] = true;
                MatFin[j][i] = true;
            }else{
                MatFin[i][j] = false;
                MatFin[j][i] = false;
            }
        }
    }
    for(i=0; i<n; i++)
        MatFin[i][i] = true;
    ok1 = true;
    while(ok1 == true){
        ok2 = true;
        for(i=1; i<n; i++){
            for(j=0; j<i; j++){
                if(MatFin[i][j] == true){
                    for(l=0; l<automat.GetV(); l++){
                        if(automat.GetV(l) != '$'){
                            p = Mat[i][l];
                            u = Mat[j][l];
                            if(p == -1 && u != -1){
                                for(k=0; k<automat.GetF(); k++)
                                    if(u == automat.GetFvector(k)){
                                        MatFin[i][j] = false;
                                        MatFin[j][i] = false;
                                        ok2 = false;
                                    }
                            }
                            if(u == -1 && p != -1){
                                for(k=0; k<automat.GetF(); k++)
                                    if(p == automat.GetFvector(k)){
                                        MatFin[i][j] = false;
                                        MatFin[j][i] = false;
                                        ok2 = false;
                                        //cout<<2<<endl;
                                    }
                                }
                            if(p != -1 && u != -1)
                                if(MatFin[p][u] == false){
                                    MatFin[i][j] = false;
                                    MatFin[j][i] = false;
                                    ok2 = false;
                            }
                        }
                    }
                }
            }
        }
        if(ok2 == true)
            ok1 = false;
    }

    for(i=0; i<=n; i++){
        for(j=0; j<=n; j++){
            multimi[i][j] = -1;
        }
    }
    p = -1;
    for(i=n-1; i>=0; i--){
        if(fr[i] != 1){
            fr[i] = 1;
            p++;
            ///creez o multime pentru i
            multimi[p][0] = i;
            for(j=i-1; j>=0; j--){
                ok1 = true;
                if(fr[j] != 1){
                    if(MatFin[i][j] != true)
                        ok1 = false;
                }else ok1 = false;
                if(ok1 == true){
                    fr[j] = 1;
                    /// adaug j in aceeasi multime in care se afla si i
                    u = 1;
                    while(multimi[p][u] != -1)
                        u++;
                    while(u!=0){
                        multimi[p][u] = multimi[p][u-1];
                        u--;
                    }
                    multimi[p][0] = j;
                }
            }
        }
    }

    for(l=0; l<automat.GetV(); l++)
    {
        for(i=0;i<n;i++)
        {
            for(j=0; j<n; j++){
                if(automat.GetV(l) != '$')
                    VecMat[l][i][j] = -1;
            }
        }
    }

    for(l=0; l<n; l++){
        if(multimi[l][0] != -1){
            for(j=0; j<automat.GetV(); j++){
                if(automat.GetV(j) != '$'){
                    m = 0;
                    ok2 = true;
                    while(Mat[multimi[l][m]][j] == -1 && ok2 == true){
                        m++;
                        if(multimi[l][m] == -1)
                            ok2 = false;
                    }
                    if(ok2 == true){
                        p = -1;
                        ok1 = false;
                        while(ok1 == false){
                            u = 0;
                            p++;
                            while(multimi[p][u] != -1){
                                if(multimi[p][u] == Mat[multimi[l][m]][j])
                                    ok1 = true;
                                u++;
                            }
                        }
                        for(i=0; i<n; i++){
                            VecMat[j][l][i] = multimi[p][i];
                        }
                    }
                }
            }
        }
    }

    p = 0;
    for(i=0; i<n; i++){
        ok1 = false;
        for(j=0; j<n; j++){
            if(multimi[i][j] != -1){
                for(l=0; l<automat.GetF(); l++){
                    if(multimi[i][j] == automat.GetFvector(l)){
                        fr[p] = i;
                        p++;
                        ok1 = true;
                        break;
                    }
                }
            }
            if(ok1 == true)
                break;
        }
    }

    for(i=0; i<n; i++)
        for(j=0; j<n; j++)
            if(multimi[i][j] == 0)
                si = i;
    nr = 0;
    for(i=0; i<n; i++)
        if(multimi[i][0] != -1){
            frecv[nr] = 0;
            nr++;
        }

    for(i=0; i<automat.GetV(); i++){
        for(j=0; j<nr; j++){
            if(VecMat[i][j][0] != -1)
                frecv[j]++;
        }
    }
    ok2 = true;
    while(ok2 == true){
        ok2 = false;
        for(j=0; j<nr; j++)
            if(frecv[j] == 0 && j != si){
                ok1 = false;
                for(i=0; i<p; i++){
                    if(j == fr[i])
                        ok1 = true;
                }
                if(ok1 == false){
                    ok2 = true;
                    frecv[j] = 1;
                    u = multimi[j][0];
                    ///nodul de pe pozitia j din multimi trebuie sters
                    for(i=0; i<n; i++){
                        multimi[j][i] = -1;
                    }
                    for(i=0; i<automat.GetV(); i++){
                        for(l=0; l<nr; l++){
                            if(VecMat[i][l][0] == u){
                                for(k=0; k<n; k++){
                                    VecMat[i][l][k] = -1;
                                }
                            }
                        }
                    }
                }
            }
    }

    for(l=0; l<n; l++){
        ok1 = true;
        for(i=0; i<automat.GetV()-1; i++){
            if(automat.GetV(i) != '$' && automat.GetV(i+1) != '$')
                if(VecMat[l][i][0] != multimi[l][0] && VecMat[l][i][0] != -1){
                    ok1 = false;
                }
        if(l == si)
            ok1 = false;
        }
        if(ok1 == true){
            for(i=0; i<p; i++){
                if(fr[i] == l)
                    ok1 = false;
            }
            if(ok1 == true){
                u = multimi[l][0];
                for(i=0; i<n; i++){
                    multimi[l][i] = -1;
                }
                for(i=0; i<automat.GetV(); i++){
                    for(j=0; j<nr; j++){
                        if(VecMat[i][j][0] == u){
                            for(k=0; k<n; k++){
                                VecMat[i][j][k] = -1;
                            }
                        }
                    }
                }
            }
        }
    }

    nr = 0;
    u = 0;
    ok1 = true;
    for(i=0; i<n; i++){
        frecv[i] = 0;
        if(ok1 == true)
            if(multimi[i][0] != -1){
                ok1 = false;
                nr = i;
            }
        if(multimi[i][0] != -1)
            u = i;
    }
    for(i=nr; i<=u; i++){
        if(multimi[i][0] != -1){
            ok1 = true;
            for(l=0; l<n; l++)
                for(j=0; j<automat.GetV(); j++){
                    if(VecMat[j][l][0] == multimi[i][0]){
                        ok1 = false;
                    }
                }
            if(i == si)
                ok1 = false;
            if(ok1 == true){
                for(j=0; j<p; j++)
                    if(fr[j] == i){
                        for(l=j; l<p; l++)
                            fr[l] = fr[l+1];
                        p--;
                    }

                for(j=0; j<n; j++)
                    multimi[i][j] = -1;
            }
        }
    }
    k = 0;
    for(i=nr; i<=u; i++)
        if(multimi[i][0] != -1)
            k++;
    /*  ///afisare pentru ajutor
    for(l=0; l<automat.GetV(); l++)
    {
        if(automat.GetV(l) != '$'){
            for(i=0;i<n;i++)
            {
                for(j=0; j<n; j++){
                        cout<<VecMat[l][i][j]<<" ";
                }
                cout<<endl;
            }
            cout<<endl;
        }
    }
    for(i=0; i<n; i++){
        cout<<i<<": ";
        for(j=0; j<n; j++){
            cout<<multimi[i][j]<<" ";
        }
        cout<<endl;
    }
    */
    automat.SetQ(k);
    automat.SetF(p);
    for(i=0; i<p; i++){
        fr[i] = multimi[fr[i]][0];
    }
    for(i=0; i<p; i++)
        automat.SetFvector(fr, i);
    k=0;

    for(i=0; i<automat.GetV(); i++){
        if(automat.GetV(i) != '$'){
            for(j=0; j<n; j++){
                if(multimi[j][0] != -1 && VecMat[i][j][0] != -1){
                    automat.SetTriplete(multimi[j][0], automat.GetV(i), VecMat[i][j][0], k);
                    k++;
                }
            }
        }

    }
    automat.SetG(k);
    cout<<"Numarul Starilor a fost inlocuit cu prima cifra a multimilor formate din noduri echivalente."<<endl;

    return automat;
}

ostream &operator << (ostream &out, NFA &auxo)
{
    out<<"Numarul total de stari este: "<<auxo.GetQ();
    out<<endl;
    out<<"Numarul de stari finale este: "<<auxo.GetF();
    out<<endl<<"Starile finale sunt: ";
    for(int i=0; i<auxo.GetF(); i++)
        out<<auxo.GetFvector(i)<<" ";
    out<<endl;
    out<<"Numarul de elemente al alfabetului este: "<<auxo.GetV();
    out<<endl<<"Alfabetul este: ";
    for(int i=0; i<auxo.GetV(); i++)
        out<<auxo.GetV(i)<<" ";
    out<<endl;
    out<<"Numarul de tranzitii este: "<<auxo.GetG();
    out<<endl<<"Tranzitiile sunt:\n";
    for(int i=0; i<auxo.GetG(); i++)
        out<<auxo.GetGa(i)<<" "<<auxo.GetGc(i)<<" "<<auxo.GetGb(i)<<endl;
    return out;
}




int main()
{
    NFA automat("afn.in");
    automat = lnfa_to_nfa(automat);
    cout<<automat;

    cout<<endl;

    automat = nfa_to_dfa(automat);
    cout<<automat;

    cout<<endl;

    automat = dfa_to_dfamin(automat);
    cout<<automat;

    return 0;
}
