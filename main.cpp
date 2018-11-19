#include <iostream>
#include <vector>
#include <ctime>
#include <algorithm>

using namespace std;

int MutationRate = 25; // We are going to use this for mutation of genes 1000 = 100% 25 = .25% mutation rate
int CrossRate = 500;
int n , d ;
double fsum = 0;
pair<double,double> point[500];

struct Member {
    double DNA[5];
    double Fitness;
};

const int PopSize = 2000;
struct Population {
    vector<Member> Members = vector<Member>(PopSize);
};

void calcFitness(Member &member){
    member.Fitness = 0;
    for(int i=0 ; i<n ; i++){
        double yCalc = 0;
        for(int j=0 ; j<=d ; j++)
            yCalc += member.DNA[j] * pow(point[i].first,j);
        member.Fitness += (yCalc - point[i].second)*(yCalc - point[i].second);
    }
    member.Fitness /= (double)n;
}

void crossOver(Member parent1 , Member parent2 , Population &nextPop , int j){
    int randPos = rand() % d;
    for(int i = randPos+1 ; i<=d ; i++)
        swap(parent1.DNA[i],parent2.DNA[i]);
    nextPop.Members[j] = parent1;
    nextPop.Members[j+1] = parent2;
}

double randDouble(double l , double r){
    return l + (double) rand() / ((double)RAND_MAX/(r-l));
}

void Mutation(double &val , int gen , int mGen){
    double left = val - (-10.0);
    double right = (10.0) - val;
    int r1 = rand()%1000;
    double r2 = randDouble(0,1);
    double b = randDouble(0.5,5.0);
    double y = r1 < 500 ?left:right;
    double deltaMutation = y * (1 - pow(r2 , pow((1 - gen/mGen), b) ) );

    if(r1 < 500)val = val - deltaMutation;
    else val = val + deltaMutation;
}

Member rouletteWheel(Population Pop){
    double cumulative = 0;
    for(int i=0 ; i< PopSize ; i++) {
        double randSelect = randDouble(0,fsum);
        for (auto &Member : Pop.Members) {
            cumulative += Member.Fitness;
            if (cumulative > randSelect)
                return Member;
        }
    }
}


int main() {

    //freopen("input_example.txt", "r", stdin);

    srand(time(nullptr));

    int t;
    cin>>t;
    while(t--) {
        cin >> n >> d;
        for (int i = 0; i < n; i++)
            cin >> point[i].first >> point[i].second;

        //create a population and initialize it with random DNA, also set the fitness to 0
        Population Pop , nextPop , tempPop;
        for (auto &Member : Pop.Members)
            for (int j = 0; j <= d; j++)
                Member.DNA[j] = randDouble(-10,10);

        int Generation = 0 , MaxGen = 500;
        while (Generation++ < MaxGen ) {

            fsum = 0;
            for(auto &Member : Pop.Members){
                calcFitness(Member);
                fsum += Member.Fitness;
            }

            double tempFsum = 0;
            for(auto &Member : Pop.Members) {
                Member.Fitness = fsum - Member.Fitness;
                tempFsum += Member.Fitness;
            }
            fsum = tempFsum;

            sort(Pop.Members.begin(), Pop.Members.end(),
                 [](Member const &a, Member &b) { return a.Fitness > b.Fitness; });

            //Selection & CrossOver
            for(int i=0 ; i<PopSize ; i += 2){
                Member parent1 = rouletteWheel(Pop);
                Member parent2 = rouletteWheel(Pop);
                if(rand()%1000 < CrossRate)
                    crossOver(parent1,parent2,nextPop,i);
                else{
                    nextPop.Members[i] = parent1;
                    nextPop.Members[i+1] = parent2;
                }
            }

            //Mutation
            for(auto &Member : nextPop.Members){
                for(int i=0 ; i <= d; i++)
                    if(rand()%1000 < MutationRate)
                        Mutation(Member.DNA[i] , Generation , MaxGen);
                calcFitness(Member);
            }
            for(auto &Member : Pop.Members)
                calcFitness(Member);

            //now lets sort the population by fitness, from Lowest
            sort(nextPop.Members.begin(), nextPop.Members.end(),
                 [](Member const &a, Member &b) { return a.Fitness < b.Fitness; });

            sort(Pop.Members.begin(), Pop.Members.end(),
                 [](Member const &a, Member &b) { return a.Fitness < b.Fitness; });

            for(int i=PopSize/2 ; i<PopSize ; i++){
                Pop.Members[i] = nextPop.Members[i-(PopSize/2)];
            }

            sort(Pop.Members.begin(), Pop.Members.end(),
                 [](Member const &a, Member &b) { return a.Fitness < b.Fitness; });

            //lets print some stuff
            cout << "Generation : " << Generation << "\tLowest Fitness : " << Pop.Members[0].Fitness << " With Sequence : ";
            for(int i=0 ; i<=d ; i++)cout<<Pop.Members[0].DNA[i]<<" \n"[i==d];
        }
        cout<<"Generation evolved!"<<endl<<endl;
    }

    return 0;
}
