#include <iostream>
#include <string>
#include <fstream>
#include <unistd.h>
#include <ctype.h>
#include <set>
#include <map>
#include <list>
#include <vector>
#include <cmath>
#include "net.h"
#include "cell.h"
#include <ctime>

using namespace std;

typedef vector <int> vi;

vector <Net*> vn, *bestvn;
vector <Cell*> vc, *bestvc;
vector <int> cellstack;
int k = 0, bestk;
bool *bestset;
vector <int> bestA, bestB;
map <string, int> mc, mn; 
int ccnt = 0, ncnt = 0;

int cutSize = 0;
ifstream ifc, ifn;
ofstream of;
double error; 
int tcsz = 0, acsz = 0, bcsz = 0, cs = 0;
int accnt = 0, bccnt = 0, afccnt = 0, bfccnt = 0, accg = 0, bestg = 0;

int bestacnnt = 0, bestbcnnt = 0, bestacsz = 0, bestbcsz = 0;
int Pmax = 0;
double tstart, tend;
map <int, Node*> blist[2];

void parseCells(istream & in){
    string str;
    int size;
    while (in >> str >> size){
        mc[str] = ccnt;
        if (acsz <= bcsz){
            Cell *c = new Cell(str, size, 0, ccnt);
            vc.push_back(c);
            acsz += size;
            accnt++;
        }
        else {
            Cell *c = new Cell(str, size, 1, ccnt);
            vc.push_back(c);
            bcsz += size;
            bccnt++;
        }
        tcsz += size; 
        ccnt++;
    }
}

void calAB(){
    for (int i = 0; i < ncnt; i++){
        vi & vl = vn[i]->cellList;
        vn[i]->B = 0;
        vn[i]->A = 0;
        for (int j = 0; j < vl.size(); j++){
            Cell * cell = vc[vl[j]];
            if (cell->set) vn[i]->B++;
            else vn[i]->A++;
        }
    }
}

void parseNets(istream & in){
    string str, tmp;
    while (in >> tmp){ // NET
        in >> str;  // nxxx
        mn[str] = ncnt;
        in >> tmp;  // {
        Net *n = new Net(str);
        vn.push_back(n);
        while (in >> tmp && tmp[0] != '}'){
            vi & l = vc[mc[tmp]]->netList;
            if (!l.size() || l[l.size()-1] != ncnt) {
                l.push_back(ncnt);
                vc[mc[tmp]]->pins++;
                vn[ncnt]->cellList.push_back(mc[tmp]);
                if (vc[mc[tmp]]->set) vn[ncnt]->B++;
                else vn[ncnt]->A++;
            }
        }
        ncnt++;
    }
}

void countCutSize(){
    cs = 0;
    for (int i = 0; i < ncnt; i++)
        if (vn[i]->A && vn[i]->B) cs++;
}


void test(){
    for (int i = 0; i < ccnt; i++){
        cout << vc[i]->name << " s"
	     << vc[i]->size << " g"
             << vc[i]->gain << " ab"
	     << vc[i]->set  << " p"
	     << vc[i]->pins << " l"
	     << vc[i]->lock << ' ';
	for (int j = 0; j < vc[i]->netList.size(); j++){
	    int id = vc[i]->netList[j];
	    cout << vn[id]->name << ' ';
        }
	cout << endl;
    }
    cout << "...\n";
    cout << "# of cells = " << ccnt << endl;
    cout << "Total Cell Size = " << tcsz << endl;
    cout << "Set A Size = " << acsz << endl;
    cout << "Set B Size = " << bcsz << endl;
    cout << "|A - B| = " << abs(acsz-bcsz) << endl;
    cout << "...\n";
    for (int i = 0; i < ncnt; i++){
        cout << vn[i]->name << ' ';
        for (int j = 0; j < vn[i]->cellList.size(); j++){
            int id = vn[i]->cellList[j];
            cout << vc[id]->name << ' ';
        }
        cout << endl;
    }
    cout << "...\n";
    return;
}

void countPmax(){
    for (int i = 0; i < ccnt; i++)
        if (vc[i]->pins > Pmax) Pmax = vc[i]->pins;
}

void countError(){
    error = (double) tcsz/10;
}

void outputFile(ostream & out){
    out << "cut_size " << cs << endl;
    out << "A " << accnt << endl;
    for (int i = 0; i < ccnt; i++)
        if (!vc[i]->set)
            out << vc[i]->name << endl;
    out << "B " << bccnt << endl;
    for (int i = 0; i < ccnt; i++)
        if (vc[i]->set)
            out << vc[i]->name << endl;
}


void parseInput(int argc, char ** argv){

    char opt = 0;
    while ((opt = getopt(argc, argv, "c:n:o:h?")) != -1){
        switch (opt){
            case 'c':
                ifc.open(optarg, ios::in);
                if (!ifc.is_open())
                    cout << "Cannot open the cells file at [-" << opt << ' ' << optarg << ']' << endl;		
                break;
            case 'n':
                ifn.open(optarg, ios::in);
                if (!ifn.is_open())
                    cout << "Cannot open the nets file at [-" << opt << ' ' << optarg << ']' << endl;		
                break;
            case 'o':
                of.open(optarg, ios::out);
                if (!of.is_open()) {
                    cout << "Cannot open the output file at [-" << opt << ' ' << optarg << ']' << endl;		
                }
                break;
            case 'h' :
			case '?' :
            default:
                cerr << "Usage: " << argv[0] << " -c <cells file name> -n <nets file name> -o <output file name>\n";
                exit(EXIT_FAILURE);
        }
    }

}


void traverse(){
    
    for (int k = 0; k < 2; k++){
        cout << "---- " << ((!k) ? "A" : "B") << " ----\n";
        for (int i = Pmax ; i >= -Pmax; i--){
            cout << '[' << i << ']' << ' ';
            Node *trav = blist[k][i]->next;
            while (trav != NULL){
                cout << vc[trav->id]->name << "->";
                trav = trav->next;
            }
            cout << endl;
        }
    }


}

void remove(Cell * c){
    Node *p = c->to;
    p->prev->next = p->next;
    if (p->next != NULL) p->next->prev = p->prev;
}


void insert_front(Cell * c){
    int gain = c->gain;
    bool set = c->set;
    Node *p = c->to;
    p->prev = blist[set][gain];
    p->next = blist[set][gain]->next;
    blist[set][gain]->next = p;
    if (p->next != NULL) p->next->prev = p;
}

void move(Cell * c){
    remove(c);
    insert_front(c);
}

void buildBlist(){
    blist[0].clear();
    blist[1].clear();
    for (int i = -Pmax; i <= Pmax; i++) {
        if (blist[0][i] == NULL) blist[0][i] = new Node(-1);
        if (blist[1][i] == NULL) blist[1][i] = new Node(-1);
    }
    for (int i = 0; i < ccnt; i++)
        insert_front(vc[i]);
}


Cell * findMaxGain(bool set){
    int p = Pmax;
    while (p >= -Pmax && blist[set][p]->next == NULL){p--;}
    Cell * ans = vc[blist[set][p]->next->id];
    return ans;
}


void reverse(){
    int i = cellstack.size()-1;
    for (; i > k; i--)
        //cout << cellstack[i] << " ";
        vc[cellstack[i]]->set = !vc[cellstack[i]]->set;
    
}

void store(){
    bestg = accg;
    //bestvc = &vc;
    //bestvn = &vn;
    bestacnnt = accnt;
    bestbcnnt = bccnt;
    bestacsz = acsz;
    bestbcsz = bcsz;
    //bestset.clear();
    //bestA.clear();
    //bestB.clear();
    bestk = k;
    //for (int i = 0 ; i < ccnt; i++)
    //    bestset[i] = vc[i]->set;
}



void restore(){
    //vc = *bestvc;
    //vn = *bestvn;
    k = bestk;
    //cout << k;
    accnt = bestacnnt;
    bccnt = bestbcnnt;
    acsz = bestacsz;
    bcsz = bestbcsz;
    //for (int i = 0 ; i < ccnt; i++)
    //    vc[i]->set = bestset[i];
    reverse();
    calAB();
    //cout << "???\n";
}


void initGain(){
    for (int i = 0; i < ccnt; i++){
        vc[i]->gain = 0;
        vc[i]->lock = 0;
    }
    
    accg = 0;
    store();
    afccnt = accnt;
    bfccnt = bccnt;
   
    
    
    for (int i = 0; i < ccnt; i++){
        for (int j = 0 ; j < vc[i]->netList.size(); j++){
            int id = vc[i]->netList[j];
            if (vc[i]->set == 0) {
                if (vn[id]->A == 1) vc[i]->gain++;
                if (vn[id]->B == 0) vc[i]->gain--;
            }
            else {
                if (vn[id]->B == 1) vc[i]->gain++;
                if (vn[id]->A == 0) vc[i]->gain--;
            }
        }
    }
    buildBlist();
}

void updateGain(Cell * c){
    accg += c->gain;

    
    c->lock = true;
    int num = c->to->id;
    cellstack.push_back(num);
    if (!c->set) {
        int szn = c->netList.size();
        for(int i = 0; i < szn; i++){
            int id = c->netList[i];
            Net * net = vn[id];
            int szc = net->cellList.size();
            if (net->B == 0){
                for (int j = 0; j < szc; j++){
                    int idc = net->cellList[j];
                    if (!vc[idc]->lock) {
                        vc[idc]->gain++;
                        move(vc[idc]);
                    }
                }
            }
            else if (net->B == 1){
                for (int j = 0; j < szc; j++){
                    int idc = net->cellList[j];
                    if (!vc[idc]->lock && vc[idc]->set) {
                        vc[idc]->gain--;
                        move(vc[idc]);
                    }
                }
            }
            net->A--;
            net->B++;
            c->set = true;
            if (net->A == 0){
                for (int j = 0; j < szc; j++){
                    int idc = net->cellList[j];
                    if (!vc[idc]->lock) {
                        vc[idc]->gain--;
                        move(vc[idc]);
                    }
                }
            }
            else if (net->A == 1){
                for (int j = 0; j < szc; j++){
                    int idc = net->cellList[j];
                    if (!vc[idc]->lock && !vc[idc]->set) {
                        vc[idc]->gain++;
                        move(vc[idc]);
                    }
                }
            }
        }
        remove(c);
        acsz -= c->size;
        bcsz += c->size;
        afccnt--;
        accnt--;
        bccnt++;
    }
    else {
        int szn = c->netList.size();
        for(int i = 0; i < szn; i++){
            int id = c->netList[i];
            Net * net = vn[id];
            int szc = net->cellList.size();
            if (net->A == 0){
                for (int j = 0; j < szc; j++){
                    int idc = net->cellList[j];
                    if (!vc[idc]->lock) {
                        vc[idc]->gain++;
                        move(vc[idc]);
                    }
                }
            }
            else if (net->A == 1){
                for (int j = 0; j < szc; j++){
                    int idc = net->cellList[j];
                    if (!vc[idc]->lock && !vc[idc]->set) {
                        vc[idc]->gain--;
                        move(vc[idc]);
                    }
                }
            }
            net->B--;
            net->A++;
            c->set = false;
            if (net->B == 0){
                for (int j = 0; j < szc; j++){
                    int idc = net->cellList[j];
                    if (!vc[idc]->lock) {
                        vc[idc]->gain--;
                        move(vc[idc]);
                    }
                }
            }
            else if (net->B == 1){
                for (int j = 0; j < szc; j++){
                    int idc = net->cellList[j];
                    if (!vc[idc]->lock && vc[idc]->set) {
                        vc[idc]->gain++;
                        move(vc[idc]);
                    }
                }
            }
        }
        remove(c);
        bcsz -= c->size;
        acsz += c->size;
        bfccnt--;
        bccnt--;
        accnt++;
    }
    if (accg > bestg)
        store();
    return;
}

int pass = 0;

void FMAlgorithm(){
    bool flag = false;
    initGain();
    int count = 0;
    k = 0;
    bestk = 0;
    cellstack.clear();
    while (!flag && count++ < ccnt){
        if (!bfccnt){
            Cell * a = findMaxGain(0);
            if (abs(acsz-bcsz-2*a->size) < error) updateGain(a);
            else flag = true;
        }
        else if (!afccnt){
            Cell * b = findMaxGain(1);
            if (abs(bcsz-acsz-2*b->size) < error) updateGain(b);
            else flag = true;
        }
        else {
            Cell * a = findMaxGain(0), * b = findMaxGain(1);
            if (a->gain >= b->gain) {
                if (abs(acsz-bcsz-2*a->size) < error) updateGain(a);
                else if (abs(bcsz-acsz-2*b->size) < error) updateGain(b);
                else flag = true;
            }
            else {
                if (abs(bcsz-acsz-2*b->size) < error) updateGain(b);
                else if (abs(acsz-bcsz-2*a->size) < error) updateGain(a);
                else flag = true;
            }
        }
        k++;
    }
    
    if (bestg > 0 ) {
        pass++;
        
        restore();
        cout << "Pass " << pass << endl;
        cout << "Best Partial Sum of Gains: " << bestg << endl;
        cout << "Total Sum of Gains (Should be 0): " << accg << endl;
        cout << endl;
        FMAlgorithm();
        
    }
    else { bestk = -1; k = -1; return;}
}

void adjust(){
    if (abs(acsz-acsz) < error) return;
    else {
        cout << "...Need balancing.\n";
        int i;
        for (i = 0; i < ccnt && abs(acsz-bcsz) >= error; i++){
            Cell * c = vc[i];
            if (acsz > bcsz && !c->set){
                acsz -= c->size;
                bcsz += c->size;
                c->set = true;
            }
            else if (acsz < bcsz && c->set){
                acsz += c->size;
                bcsz -= c->size;
                c->set = false;
            }
        }
        if (i == ccnt && abs(acsz-bcsz) >= error) {
            cerr << "(ERROR)...This testcase can never be balanced!\n";
            //exit(EXIT_FAILURE);
        }
    }
}


int main(int argc, char *argv[]){
    
    ios_base::sync_with_stdio(false);
	
    parseInput(argc, argv);
    if (ifc.is_open()) parseCells(ifc);
    else parseCells(cin);
    ifc.close();

    if (ifn.is_open()) parseNets(ifn);
    else parseNets(cin);
    ifn.close();
    //bestset = new bool[ccnt]();
    countCutSize();
    cout << "Initial Cut Size = " << cs << endl;
    cout << endl;
    countError();
    countPmax();
    adjust();
    
    bestacsz = acsz;
    bestbcsz = bcsz;
    tstart = clock();
    FMAlgorithm();
	tend = clock();
    restore();
    countCutSize();
    cout << "Final Cut Size = " << cs << endl;
    if (of.is_open()) outputFile(of);
    else outputFile(cout);
    of.close();
    cout << endl;
    cout << "FM Algorithm Run Time: " << (double)(tend-tstart)/CLOCKS_PER_SEC << " sec\n";
    cout << "Total Run Time: " << (double)clock()/CLOCKS_PER_SEC << " sec\n";
    
}
