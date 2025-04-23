#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <bits/stdc++.h>
using namespace std;

struct Jakobian;
struct ElemUniv;


struct node {
   friend ostream& operator<<(ostream& str, const node& n) {
       str << "x: " << n.x << " y: " << n.y << " BC: " << (n.BC ? "true" : "false");
       return str;
    }
    float x{}, y{};
    bool BC{};


    node(const float x, const float y) : x(x), y(y) {
        this -> x = x;
        this -> y = y;
        BC = false;
    }

    node() = default;

    void set_bc_condition(bool condition) {
        BC = condition;
    }

};

vector<string> split(const string &str, const char del) {
    vector<string> tokens;
    stringstream ss(str);
    string word;
    while (!ss.eof()) {
        getline(ss, word, del);
        tokens.push_back(word);
    }
    return tokens;
}


    struct GlobalData {
        int SimulationTime;
        int SimulationStepTime{};
        int Conductivity{};
        int Alfa{};
        int Tot{};
        int InitialTemp{};
        int Density{};
        int SpecificHeat{};
        int nN{};
        int nE{};

        explicit GlobalData(const string &fileName_) {
            int variables[10];
            ifstream MyReadFile("..\\" + fileName_);
            string MyText;
            int i = 0;
            while (getline(MyReadFile, MyText) && i < 10) {
                vector<string> tokens = split(MyText, ' ');
                variables[i] = stoi(tokens.at(tokens.size() - 1));
                i++;
            }
            this->SimulationTime = variables[0];
            this->SimulationStepTime = variables[1];
            this->Conductivity = variables[2];
            this->Alfa = variables[3];
            this->Tot = variables[4];
            this->InitialTemp = variables[5];
            this->Density = variables[6];
            this->SpecificHeat = variables[7];
            this->nN = variables[8];
            this->nE = variables[9];
        }

        void PrintData() const {
            cout << "Simulation Time: " << SimulationTime << endl;
            cout << "Simulation Step Time: " << SimulationStepTime << endl;
            cout << "Conductivity: " << Conductivity << endl;
            cout << "Alfa: " << Alfa << endl;
            cout << "Tot: " << Tot << endl;
            cout << "SpecificHeat: " << SpecificHeat << endl;
            cout << "NN: " << nN << endl;
            cout << "NE: " << nE << endl;
        }

        int getnN() const {
            return nN;
        }

        int getnE() const {
            return nE;
        }
    };

struct Surface {
    vector<vector<float>> N;
    int npc_bc{};

    Surface(vector<vector<float>> N) {
        this -> N = N;
    }

};



    struct ElemUniv {
        vector<vector<float>> dNdN;
        vector<vector<float>> dNdE;
        vector<vector<float>> N_c;
        vector <double> w;

        int npc{}, npc_bc{};

        vector <Surface> surface;


        ElemUniv(int npc, int npc_bc) {
            this->npc = npc;
            this->npc_bc = npc_bc;

            vector<vector<float>> dNdE (npc, vector<float> (4,0));
            vector<vector<float>> dNdN (npc, vector<float> (4,0));
            vector<vector<float>> N (sqrt(npc_bc), vector<float> (4));
            vector<vector<float>> N_c (npc, vector<float> (4,0));


            this->dNdN = dNdN;
            this->dNdE = dNdE;
            this->N_c = N_c;
            vector<vector <double>> pc;

            vector<vector <double>> pc_s;

            vector<float> w (npc,0);


            if (npc==4) {  //x = E, y = N

                pc = {{- 1/sqrt(3), - 1/sqrt(3)}, {1/sqrt(3), - 1/sqrt(3)}, {-1/sqrt(3), 1/sqrt(3)}, {1/sqrt(3), 1/sqrt(3)}};
                w = {1,1,1,1};
            }

            if (npc == 9) {
                pc = {{- sqrt(3.0/5.0), -sqrt(3.0/5.0)}, {0, - sqrt(3.0/5.0)}, {sqrt(3.0/5.0), - sqrt(3.0/5.0)}, {-sqrt(3.0/5.0), 0}, {0,0}, {sqrt(3.0/5.0),0}, {-sqrt(3.0/5.0), sqrt(3.0/5.0)}, {0, sqrt(3.0/5.0)}, {sqrt(3.0/5.0), sqrt(3.0/5.0)}};
                w = {5.0/9.0*5.0/9.0, 8.0/9.0*5.0/9.0, 5.0/9.0*5.0/9.0, 5.0/9.0*8.0/9.0, 8.0/9.0*8.0/9.0, 8.0/9.0*5.0/9.0, 5.0/9.0*5.0/9.0, 8.0/9.0*5.0/9.0,5.0/9.0*5.0/9.0};
            }

            if (npc == 16) {
                double a = sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0));
                double b = sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0));
                double w1 = (18.0 - sqrt(30.0)) / 36.0;
                double w2 = (18.0 + sqrt(30.0)) / 36.0;

                pc = {
                    {-a, -a}, {-a, -b}, {-a, b}, {-a, a},
                    {-b, -a}, {-b, -b}, {-b, b}, {-b, a},
                    {b, -a}, {b, -b}, {b, b}, {b, a},
                    {a, -a}, {a, -b}, {a, b}, {a, a}
                };
                w = {
                    w1 * w1, w1 * w2, w1 * w2, w1 * w1,
                    w2 * w1, w2 * w2, w2 * w2, w2 * w1,
                    w1 * w2, w2 * w2, w2 * w2, w1 * w2,
                    w1 * w1, w1 * w2, w2 * w1, w1 * w1
                };
            }

            for (int i = 0; i < npc; i++) {
                this->dNdE[i][0] = -0.25 * (1 - pc[i][1]);
                this->dNdE[i][1] = 0.25 * (1 - pc[i][1]);
                this->dNdE[i][2] = 0.25 * (1 + pc[i][1]);
                this->dNdE[i][3] = -0.25 * (1 + pc[i][1]);

                this->dNdN[i][0] = -0.25 * (1 - pc[i][0]);
                this->dNdN[i][1] = -0.25 * (1 + pc[i][0]);
                this->dNdN[i][2] = 0.25 * (1 + pc[i][0]);
                this->dNdN[i][3] = 0.25 * (1 - pc[i][0]);

                this->N_c[i][0] = 0.25 * (1 - pc[i][0])*(1 - pc[i][1]);
                this->N_c[i][1] = 0.25 * (1 + pc[i][0])*(1 - pc[i][1]);
                this->N_c[i][2] = 0.25 * (1 + pc[i][0])*(1 + pc[i][1]);
                this->N_c[i][3] = 0.25 * (1 - pc[i][0])*(1 + pc[i][1]);
            }

            // cout<<"!!!!!!!!!!!!!!!!!!N_C:"<<endl;
            // for (int i = 0; i < npc; i++) {
            //     for (int j = 0; j < 4; j++) {
            //     cout<<this->N_c[i][j]<<" ";
            //     }
            //     cout<<endl;
            // }

            // cout << "dziala" << endl;
            // cout << "dNdE: "<<endl;
            //
            // for (int i = 0; i < npc; i++) {
            //     cout << this->dNdE[i][0] << " " << this->dNdE[i][1] << " " << this->dNdE[i][2] << " " << this->dNdE[i][3] << endl;
            // }
            //
            // cout <<endl<< "dNdn: "<<endl;
            //
            // for (int i = 0; i < npc; i++) {
            //     cout << this->dNdN[i][0] << " " << this->dNdN[i][1] << " " << this->dNdN[i][2] << " " << this->dNdN[i][3] << endl;
            // }
            // cout<<endl;

            // cout << "npc_s: " << npc_bc << endl;

            if (npc_bc == 4) {
                pc_s = { {-1.0/sqrt(3.0), -1.0}, {1.0/sqrt(3.0), -1.0},
                            {1.0, -1.0/sqrt(3.0)}, {1.0, 1.0/sqrt(3.0)},
                            {1.0/sqrt(3.0), 1.0}, {-1.0/sqrt(3.0), 1.0},
                            {-1.0, 1.0/sqrt(3.0)}, {-1.0, -1.0/sqrt(3.0)}};
            }
            if (npc_bc == 9) {
                double x_pom = 3.0/5.0;
                pc_s = {
                    {-sqrt(x_pom), -1.0}, {0.0, -1.0}, {sqrt(x_pom), -1.0},
                    {1.0, -sqrt(x_pom)}, {1.0, 0.0}, {1.0, sqrt(x_pom)},
                    {sqrt(x_pom), 1.0}, {0.0, 1.0}, {-sqrt(x_pom), 1.0},
                    {-1.0, sqrt(x_pom)}, {-1.0, 0.0}, {-1.0, -sqrt(x_pom)}
                };
            }
            if (npc_bc == 16) {
                double a = sqrt(3.0/7.0 + 2.0/7.0 * sqrt(6.0/5.0));
                double b = sqrt(3.0/7.0 - 2.0/7.0 * sqrt(6.0/5.0));
                pc_s = {
                    {-a, -1.0}, {-b, -1.0}, {b, -1.0}, {a, -1.0},
                    {1.0, -a}, {1.0, -b}, {1.0, b}, {1.0, a},
                    {a, 1.0}, {b, 1.0}, {-b, 1.0}, {-a, 1.0},
                    {-1.0, a}, {-1.0, b}, {-1.0, -b}, {-1.0, -a}
                };
            }

            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < sqrt(npc_bc); j++) {

                    N[j][0] = 0.25 * (1 - pc_s[j +(i*sqrt(npc_bc))][0]) * (1 - pc_s[j+(i*sqrt(npc_bc))][1]);
                    N[j][1] = 0.25 * (1 + pc_s[j+(i*sqrt(npc_bc))][0]) * (1 - pc_s[j+(i*sqrt(npc_bc))][1]);
                    N[j][2] = 0.25 * (1 + pc_s[j+(i*sqrt(npc_bc))][0]) * (1 + pc_s[j+(i*sqrt(npc_bc))][1]);
                    N[j][3] = 0.25 * (1 - pc_s[j+(i*sqrt(npc_bc))][0]) * (1 + pc_s[j+(i*sqrt(npc_bc))][1]);
                }
                Surface s = Surface(N);
                this->surface.push_back(s);
            }



            // cout<<"SURFACE N "<<endl;
            // cout<< surface.size()<<endl;
            //
            // for (size_t i = 0; i < surface.size(); ++i) {
            //     for (size_t j = 0; j < surface[i].N.size(); ++j) {
            //         for (size_t k = 0; k < surface[i].N[j].size(); ++k) {
            //             cout << surface[i].N[j][k] << " ";
            //         }
            //         cout << endl;
            //     }
            //     cout << "------------------------" << endl;
            // }
            // cout << endl;

        }

        vector<vector<float>> get_dNdE() {return dNdE;}
        vector<vector<float>> get_dNdN() {return dNdN;}
        vector<double> get_w(){return w;}

    };

 struct Jakobian {
        vector<vector<float>> J;
        vector<vector<float>> J1;
        float detJ{};

        Jakobian (vector<float> dNdE, vector<float> dNdN, vector<vector<float>> fin_element) {
            vector<vector<float>> J(2, vector<float>(2,0));
            this->J = J;

            for (int i = 0; i < 2; i++) {
                // cout<<endl<<i<<endl;
                for (int j = 0; j < 4; j++) {
                    // cout<<endl;
                    // cout<<"dNdE: "<<endl;
                    // cout<<fin_element[i][j]<<endl;
                    // cout<<dNdE[j]<< " * "<< fin_element[i][j]<<endl;;
                    // cout<<"dNdN: "<<endl;
                    // cout<<dNdN[j]<< "[i] * "<< fin_element[i][j]<<endl;

                    this->J[0][i] += dNdE[j] *fin_element[i][j];
                    // cout<<J[0][i]<<" = "<< dNdE[j]<<" * " << fin_element[i][j]<<endl;

                    this->J[1][i] += dNdN[j] *fin_element[i][j];
                }
            }

            // cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!LALALALALALALALA"<<endl;
            // cout<<endl;
            // for (int i =0; i<2;i++) {
            //         cout<<this->J[0][i]<<" "<<this->J[1][i] <<endl;
            // }

            ///DET J
            detJ = this->J[0][0] * this->J[1][1] - this->J[0][1] * this->J[1][0];

            // cout<<endl<<"det J: "<<detJ<<endl;
            ///INVERTED JACOBIAN
            vector<vector<float>> J1(2, vector<float>(2,0));
            this->J1 = J1;

            this->J1[0][0] = this->J[1][1];
            this->J1[0][1] = -this->J[0][1];
            this->J1[1][0] = -this->J[1][0];
            this->J1[1][1] =this-> J[0][0];
        }

        void printJacobian() const {
            cout << "Jacobian:\n";
            cout << "[ " << J[0][0] << ", " << J[0][1] << " ]\n";
            cout << "[ " << J[1][0] << ", " << J[1][1] << " ]\n";
            cout << "Det J: " << detJ << "\n";

            // cout << "Inverted Jacobian:\n";
            // cout << "[ " << J1[0][0] << ", " << J1[0][1] << " ]\n";
            // cout << "[ " << J1[1][0] << ", " << J1[1][1] << " ]\n";
        }


     vector<vector<float>> get_J() {return J;};
     vector<vector<float>> get_J1() {return J1;}
     float get_detJ() {return detJ;}
    };


struct RozwiazanieUkladuRownan {
    vector<vector<float>> H_global;            //tu zmienic// czemó??// jest db
    vector<float> P_global;
    vector<vector<float>> C_global;
    vector<vector<float>> HC;
    vector<float> PC;
    vector<float> t0;



    RozwiazanieUkladuRownan(int nN, int InitialTemp) {
        H_global = vector<vector<float>>(nN, vector<float>(nN, 0));
        P_global = vector<float>(nN, 0);
        C_global = vector<vector<float>>(nN, vector<float>(nN, 0));
        HC = vector<vector<float>>(nN, vector<float>(nN, 0));
        PC = vector<float>(nN, 0);
        t0 = vector<float>(nN, InitialTemp);

    }

    void set_to(vector<float>t0) {
        this->t0 = t0;
    }

    void print_H_global() {
        cout << endl << "GLOBAL H: " << endl;
        for (int i = 0; i < H_global.size(); i++) {
            for (int j = 0; j < H_global[i].size(); j++) {
                cout << H_global[i][j] << " ";
            }
            cout << endl;
        }
    }

    void print_HC() {
        cout << endl << "HC: " << endl;
        for (int i = 0; i < HC.size(); i++) {
            for (int j = 0; j < HC[i].size(); j++) {
                cout << HC[i][j] << " ";
            }
            cout << endl;
        }
    }

    void print_P_global() {
        cout << endl << "GLOBAL P: " << endl;
        for (int i = 0; i < P_global.size(); i++) {
            cout << P_global[i] << " ";
        }
    }

    void print_t0() {
        cout << endl << "T0: " << endl;
        for (int i = 0; i < t0.size(); i++) {
            cout << t0[i] << " ";
        }
    }

    void print_PC() {
        cout << endl << "PC: " << endl;
        for (int i = 0; i < PC.size(); i++) {
            cout << PC[i] << " ";
        }
    }

    void print_C_global() {
        cout << endl << "GLOBAL C: " << endl;
        for (int i = 0; i < C_global.size(); i++) {
            for (int j = 0; j < C_global[i].size(); j++) {
                cout << C_global[i][j] << " ";
            }
            cout << endl;
        }
    }


    void clear_HC() {
        for (int i = 0; i < HC.size(); i++) {
            for (int j = 0; j < HC[i].size(); j++) {
                HC[i][j]=0;
            }
        }

    }
    void clear_PC() {
            for (int i = 0; i < PC.size(); i++) {
                PC[i] = 0;
            }
    }


    vector<float> Jordan_method(vector<vector<float>> HC, vector<float> PC) {
        int n = HC.size();

        for (int i = 0; i < n; ++i) {
            HC[i].push_back(PC[i]);
        }

        for (int i = 0; i < n; ++i) {
            double diagElement = HC[i][i];
            for (int j = 0; j <= n; ++j) {
                HC[i][j] /= diagElement;
            }

            for (int k = 0; k < n; ++k) {
                if (i != k) {
                    double factor = HC[k][i];
                    for (int j = 0; j <= n; ++j) {
                        HC[k][j] -= factor * HC[i][j];
                    }
                }
            }
        }

        vector<float> u(n);
        for (int i = 0; i < n; ++i) {
            u[i] = HC[i][n];
        }

        return u;
    }

};


struct element {
    int npc{}, npc_bc{};
    vector<Jakobian> jakobian{};
    vector<vector<float>> finiteElement;
    int ID[4]{};

    friend ostream &operator<<(ostream &str, const element &e) {
        str << e.ID[0] << " " << e.ID[1] << " " << e.ID[2] << " " << e.ID[3];
        return str;
    }

    element(const int x, const int y, const int z, const int w, int npc, int npc_s) {
        ID[0] = x;
        ID[1] = y;
        ID[2] = z;
        ID[3] = w;
        finiteElement = vector<vector<float>>(2);

        this->npc=npc;
        this->npc_bc=npc_s;
    }

    void set_J(ElemUniv elem) {

        for (int i = 0; i < npc; i++) {                 //  STINKY? juz g
            jakobian.push_back(Jakobian(elem.get_dNdE()[i], elem.get_dNdN()[i], finiteElement));
        }
    }

    void make_finiteElement (vector<node> nodes) {
        for (int i = 0; i < 4; i++) {
            // cout<<"RATATATA"<<endl;
            // cout<< nodes[this->ID[i]-1].x <<endl;
            finiteElement[0].push_back(nodes[this->ID[i]-1].x);
            finiteElement[1].push_back(nodes[this->ID[i]-1].y);
        }
    }

    void Calc_all(ElemUniv elem, int conduct, int alpha, vector<node> nodes, RozwiazanieUkladuRownan &rozw, int tot, int specific_heat, int density, int timestep) {
        vector<vector<float>> dNdX(npc, vector<float> (4,0));
        vector<vector<float>> dNdY(npc, vector<float> (4,0));
        vector<vector<float>> dNdN = elem.get_dNdE();
        vector<vector<float>> dNdE = elem.get_dNdN();

        vector<vector<float>> H_bc(4, vector<float>(4,0));
        vector<vector<float>> H(npc, vector<float>(4,0));

        vector<double> P_lok(4,0);
        vector<vector<float>> C_lok(4, vector<float>(4,0));

        vector<vector<float>> bc_edges(4, vector<float>(2,0));



        //tutaj ogarniamy tablice dNdX i dNdY, ktore sa [npc x 4]
        for (int i = 0; i < npc; i++) {
            float detj = 1/jakobian[i].get_detJ();
            cout<< "i: "<< i<<endl;
            cout<<"RRRRR"<<endl;
            for (int j = 0; j < 4; j++) {
                dNdX[i][j] = jakobian[i].J1[0][0] * detj * dNdN[i][j] +  jakobian[i].J1[0][1] * detj * dNdE[i][j];
                dNdY[i][j] = jakobian[i].J1[1][0] * detj * dNdN[i][j] +  jakobian[i].J1[1][1] * detj * dNdE[i][j];
                // cout<<dNdX[i][j]<<" = "<<jakobian[i].J1[0][0]<<" * "<< detj << " * " << dNdN[i][j]<<" + "
                // << jakobian[i].J1[0][1]<< " * "<< detj<<" * "<< dNdE[i][j]<<endl;
                // cout<<dNdX[i][j]<<" ";

            }
            // cout<<endl;
        }
        int npc_sqrt = (int)sqrt(npc);

        vector<double> w(npc_sqrt);
        vector<double> wei(npc);

        if (npc ==4) {
            wei={1.0,1.0,1.0,1.0};
        }
        else if (npc == 9) {
            wei = {5.0/9.0*5.0/9.0, 8.0/9.0*5.0/9.0, 5.0/9.0*5.0/9.0, 5.0/9.0*8.0/9.0, 8.0/9.0*8.0/9.0, 8.0/9.0*5.0/9.0, 5.0/9.0*5.0/9.0, 8.0/9.0*5.0/9.0,5.0/9.0*5.0/9.0};


        }
        else if (npc == 16) {
            double w1 = (18.0 - sqrt(30.0)) / 36.0;
            double w2 = (18.0 + sqrt(30.0)) / 36.0;

            wei = {
                w1 * w1, w1 * w2, w1 * w2, w1 * w1,
                w2 * w1, w2 * w2, w2 * w2, w2 * w1,
                w1 * w2, w2 * w2, w2 * w2, w1 * w2,
                w1 * w1, w1 * w2, w2 * w1, w1 * w1
            };
        }

        if (npc_bc ==4) {
            w= {1.0,1.0};

        }
        else if (npc_bc == 9) {
            w[0] = 5.0/9.0;
            w[1] = 8.0/9.0;
            w[2] = 5.0/9.0;

        }
        else if (npc_bc == 16) {
            w[0] = (18 - sqrt(30)) / 36;
            w[1] = (18 + sqrt(30)) / 36;
            w[2] = (18 + sqrt(30)) / 36;
            w[3] = (18 - sqrt(30)) / 36;
        }
        ////////////////////////////////////////////////////////// h_bc
        for (int i = 0; i < 4; i++) {
            float detJ = (sqrt(pow(nodes[this->ID[i]-1].x - nodes[this->ID[(i+1)%4]-1].x,2)+pow(nodes[this->ID[i]-1].y - nodes[this->ID[(i+1)%4]-1].y,2)))/2;
            // cout<<detJ<<endl;

            if(nodes[this->ID[i]-1].BC && nodes[this->ID[(i+1)%4]-1].BC) {
                for (int j = 0; j < 4; j++) {
                    // cout<<"jestem tutaj"<<endl;
                    for (int k = 0; k < 4; k++) {

                        // if (i >= elem.surface.size()) {
                        //     cerr << "bledny indeks surface!" << endl;
                        //     break;
                        // }

                        if(npc_bc == 4) {
                            // cout<<"jestem w hbc"<<endl;
                            // cout<< elem.surface.size()<<endl;
                            // cout<< elem.surface[i].N[0][j]<<endl;

                            H_bc[j][k] +=  detJ * alpha * ((w[0] * elem.surface[i].N[0][j]*elem.surface[i].N[0][k])+(w[1]* elem.surface[i].N[1][j]*elem.surface[i].N[1][k]));
                            // cout << "(("<< w[0] << " * " << alpha << " * "<< elem.surface[i].N[0][k] << " * " << elem.surface[i].N[0][j] << ") + ("
                            // << w[1] << " * " << alpha << " * "<< elem.surface[i].N[1][k] << " * " << elem.surface[i].N[1][j] << ")) * "<< detJ<< endl;
                        }
                        else if (npc_bc == 9) {
                            H_bc[j][k] +=  detJ * alpha * ((w[0] * elem.surface[i].N[0][j]*elem.surface[i].N[0][k])+(w[1]* elem.surface[i].N[1][j]*elem.surface[i].N[1][k])+(w[2]* elem.surface[i].N[2][j]*elem.surface[i].N[2][k]));

                        }
                        else if (npc_bc == 16) {
                            H_bc[j][k] +=  detJ * alpha *((w[0] * elem.surface[i].N[0][j]*elem.surface[i].N[0][k])+(w[1]* elem.surface[i].N[1][j]*elem.surface[i].N[1][k])+(w[2]* elem.surface[i].N[2][j]*elem.surface[i].N[2][k])+(w[3]* elem.surface[i].N[3][j]*elem.surface[i].N[3][k]));
                        }
                    }
                    // cout<< H_bc[i][j]<<" ";


                    if(npc_bc == 4) {
                        // cout<<"jestem w npc_bc"<<endl;
                        // cout<< elem.surface.size()<<endl;
                        // cout<< elem.surface[i].N[0][j]<<endl;
                        // std::cout << typeid(detJ * (alpha * (w[0] * elem.surface[i].N[0][j] * tot) + (w[1] * elem.surface[i].N[1][j] * tot))).name() << std::endl;

                        P_lok[j] +=  detJ * alpha * ((w[0] * elem.surface[i].N[0][j]*tot)+(w[1]* elem.surface[i].N[1][j]*tot));

                    }
                    else if (npc_bc == 9) {
                        P_lok[i] +=  detJ * alpha * ((w[0] * elem.surface[i].N[0][j]*tot)+(w[1]* elem.surface[i].N[1][j]*tot)+(w[2]* elem.surface[i].N[2][j]*tot));

                    }
                    else if (npc_bc == 16) {
                        P_lok[i] +=  detJ * alpha * ((w[0] * elem.surface[i].N[0][j]*tot)+(w[1]* elem.surface[i].N[1][j]*tot)+(w[2]* elem.surface[i].N[2][j]*tot)+(w[3]* elem.surface[i].N[3][j]*tot));
                    }
                }
                cout<<endl;

            }
        }

        //
        // cout<<"P lok"<<endl;
        // for (int i = 0; i < 4; i++) {
        //     cout<<P_lok[i]<<" ";
        // }
        // cout<<endl;



        vector<vector<float>> temp_x(4, vector<float> (4,0));
        vector<vector<float>> temp_y(4, vector<float> (4,0));

        // cout<<"HHHHHHHHHHHHHHHHH"<<endl;
        for (int i = 0; i < npc; i++) {
            // cout<<sth<< " = " << conduct<< " * " <<jakobian[i].get_detJ()<<endl;

            for (int j = 0; j < 4; j++) {
                // cout<<endl;

                for (int k = 0; k < 4; k++) {

                    // cout<<"dndx "<< dNdX[i][j] << " * "<<dNdX[i][k]<<endl<<endl;
                    temp_x[j][k] = dNdX[i][j]*dNdX[i][k];
                    temp_y[j][k] = dNdY[i][j]*dNdY[i][k];

                    // cout << temp_x[j][k] << " ";
                    // cout << temp_y[j][k] << endl;

                    H[j][k] +=  conduct * jakobian[i].get_detJ() * (temp_x[j][k] + temp_y[j][k]) * wei[i];
                    C_lok[j][k] += specific_heat*density*(elem.N_c[i][j]*elem.N_c[i][k])*jakobian[i].get_detJ()*wei[i];

                    // cout<<"AAAAAAAAAAAAAAAAAAAAAA"<<endl;
                    // cout<< H[j][k]<< " = " << conduct << " * ("<< temp_x[j][k] << " + " << temp_y[j][k]<<") * "<< jakobian[i].get_detJ()<<endl;
                    // cout<<(sth * (temp_x[j][k] + temp_y[j][k]))*wei[i]<<" ";
                    // cout<<H[j][k]<<" ";
                }
                cout<<endl;
            }
            cout<<endl;
        }

        cout<<endl<<"H_LOC"<<endl;

        for (int i = 0; i < 4; i++) {
            for (int k = 0; k < 4; k++) {
                cout<<H[i][k]<<" ";
            }
            cout<<endl;
        }
        // cout<<"H_LOC"<<endl<<endl;
        //
        // cout<<endl<<"HBC HBC"<<endl;
        //
        // for (int i = 0; i < 4; i++) {
        //     for (int k = 0; k < 4; k++) {
        //         cout<<H_bc[i][k]<<" ";
        //     }
        //     cout<<endl;
        // }
        // cout<<"HBC HBC"<<endl<<endl;
        //
        //
        // cout<<"C LOK"<<endl<<endl;
        // for (int i = 0; i < 4; i++) {
        //     for (int k = 0; k < 4; k++) {
        //         cout<<C_lok[i][k]<<" ";
        //     }
        //     cout<<endl;
        // }
        // cout<<"C LOK"<<endl<<endl;
        //
        //
        // cout<<endl;


        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {

                // H[j][i] += H_bc[j][i] + (H_pc[j][i] * w[w1_idx] * w[w2_idx]); //TO TAK NIE MOZE BYC, BO H_BC JEST 4x4, WIEC ZMIEN ZEBY BYLO W INNYM MIEJSCU
                // H[j][i] +=(H_pc[j][i] * w[w1_idx] * w[w2_idx]); //TO TAK NIE MOZE BYC, BO H_BC JEST 4x4, WIEC ZMIEN ZEBY BYLO W INNYM MIEJSCU

                H[i][j] += H_bc[i][j];
                // cout<< H[i][j]<<" ";

                rozw.C_global[ID[i]-1][ID[j]-1] += C_lok[i][j];
                rozw.H_global[ID[i]-1][ID[j]-1] += H[i][j]; ;
                rozw.HC[ID[i]-1][ID[j]-1] += H[i][j]+ (C_lok[i][j]/timestep);
                rozw.PC[ID[i]-1] += (C_lok[i][j]/timestep)*rozw.t0[ID[j]-1];

                // cout<<"TU HYC HYC "<<H_pc[k][i][j]<<" * "<< w[w1_idx]<<" *"<< w[w2_idx]<<endl
                // cout<<endl<<"kurka rurka"<<endl;
            }
            rozw.P_global[ID[i]-1] += P_lok[i];
            rozw.PC[ID[i]-1] += P_lok[i];

            cout<<endl;
        }
    }

    void print_finiteElement() {
        cout << "FINITE ELEMENT" << endl;
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 4; j++) {
                cout << finiteElement[i][j] << " ";
            }
            cout << endl;
        }
    }


    void print_J() {

        for (int i = 0; i < npc; i++) {
            this->jakobian[i].printJacobian();
        }
    }
};


struct grid {
        int nN{};
        int nE{};
        int npc{}, npc_s{};
        int Tot;
        int Conductivity;
        int Alpha;
        int Density;
        int SpecificHeat;
        int SimulationTimeStep;
        int SimulationTime;
        vector<node> nodes;
        vector<element> elements;
        vector<int> BC;


        grid(const string& fileName_, const int nN, const int nE, int Conductivity, int Alpha, int npc, int npc_s, int Tot, int SpecificHeat, int Density, int SimulationTimeStep, int SimulationTime) {
            this->npc = npc;
            this->npc_s = npc_s;
            this->nN = nN;
            this->nE = nE;
            this->Tot = Tot;
            this->Conductivity = Conductivity;
            this->Conductivity = Conductivity;
            this->Alpha = Alpha;
            this->SpecificHeat = SpecificHeat;
            this->Density = Density;
            this->SimulationTimeStep = SimulationTimeStep;
            this->SimulationTime = SimulationTime;
            ifstream MyReadFile("..\\" + fileName_);
            string MyText;
            int i = -1;
            while (getline(MyReadFile, MyText)) {
                i++;
                if (i < 10) {
                    continue;
                }
                if (i > 10 && i < 11 + nN) {
                    vector<string> tokens = split(MyText, ',');
                    nodes.push_back(node(stof(tokens[1]), stof(tokens[2])));
                }
                if (i > 11 + nN && i < 12 + nN + nE) {
                    vector<string> tokens = split(MyText, ',');
                    elements.push_back(element(stoi(tokens[1]), stoi(tokens[2]), stoi(tokens[3]), stoi(tokens[4]), this->npc, this->npc_s));
                }
                if (i == 12 + nN + nE + 1) {
                    vector<string> tokens = split(MyText, ',');
                    for (int i = 0; i < tokens.size(); i++) {
                        int nodeId = stoi(tokens[i]);
                        BC.push_back(stoi(tokens[i]));
                        if (nodeId > 0 && nodeId <= nodes.size()) {
                            nodes[nodeId - 1].set_bc_condition(true);
                        }
                    }
                }
            }

        }

        void setSolution(ElemUniv eu, RozwiazanieUkladuRownan &rozw) {
            for (int i = 0; i <this->SimulationTime ; i+=this->SimulationTimeStep) { //this->SimulationTime
                rozw.clear_HC();
                rozw.clear_PC();

                for (int i = 0; i <nE; i++) { //nE
                    elements[i].make_finiteElement(this->nodes);
                    // elements[i].print_finiteElement();
                    elements[i].set_J(eu);
                    // elements[i].print_J();

                    elements[i].Calc_all(eu, this->Conductivity, this->Alpha, nodes, rozw, Tot, this->SpecificHeat, this->Density, this->SimulationTimeStep);
                }
                //
                rozw.print_H_global();
                rozw.print_HC();
                cout << endl;
                rozw.print_P_global();
                rozw.print_PC();
                cout << endl;
                rozw.print_C_global();

                cout<<endl;
                cout<<"############ SOLUTION " << i+this->SimulationTimeStep<< " #############"<<endl;
                vector<float> solution = rozw.Jordan_method(rozw.HC, rozw.PC);
                rozw.set_to(solution);

                float min = *min_element(solution.begin(), solution.end());
                float max = *max_element(solution.begin(), solution.end());

                cout << "min: " << min <<" ";
                cout << "max: " << max << endl;

                for (int i = 0; i<solution.size(); i++) {
                    cout << solution[i]<<" ";
                }
            }
        }

        void printNodes() {
            for (int i = 0; i < nN; i++) {
                cout << nodes.size();
                cout << nodes[i] << endl;
            }
        }

        void printElements() {
            for (int i = 0; i < nE; i++) {
                cout << elements[i] << endl;
            }
        }
    };

int main() {
    int npc = 16;
    int npc_bc = 16;
    // GlobalData data_1("Test1_4_4.txt");
    GlobalData data_1("Test2_4_4_MixGrid.txt");
    // GlobalData data_1("Test3_31_31_kwadrat.txt");

    ElemUniv xd(npc, npc_bc);
    data_1.PrintData();
    RozwiazanieUkladuRownan rozw(data_1.nN, data_1.InitialTemp);
    // grid myGrid("Test3_31_31_kwadrat.txt", data_1.getnN(), data_1.getnE(), data_1.Conductivity, data_1.Alfa, npc, npc_bc, data_1.Tot,data_1.SpecificHeat,data_1.Density, data_1.SimulationStepTime, data_1.SimulationTime);
    grid myGrid("Test2_4_4_MixGrid.txt", data_1.getnN(), data_1.getnE(), data_1.Conductivity, data_1.Alfa, npc, npc_bc, data_1.Tot,data_1.SpecificHeat,data_1.Density, data_1.SimulationStepTime, data_1.SimulationTime);
    // grid myGrid("Test1_4_4.txt", data_1.getnN(), data_1.getnE(), data_1.Conductivity, data_1.Alfa, npc, npc_bc, data_1.Tot, data_1.SpecificHeat, data_1.Density, data_1.SimulationStepTime, data_1.SimulationTime);
    myGrid.printNodes();
    string myText;
    myGrid.printElements();
    myGrid.setSolution(xd, rozw);

    return 0;
}






