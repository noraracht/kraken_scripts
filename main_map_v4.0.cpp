// main.cpp -- defining, prototyping, and calling a function

// features:
// change signatures from 64 to 32 bit int
// replace map with arrays

// to compile:
// g++ main_map_v4.0_cp.cpp -std=c++11 -o main_map

// to run:
//./main_map -i G000016385_neighbors_copy2.fa -d 3  -o G000016385_neighbors_copy2.fa-map
//./main_map -i G000399765_fna_neighborhood_k32C.fa -d 3 -o G000399765_fna_neighborhood_k32C-map
//./main_map -i G000307305_dump.fa -d 3 -o G000307305_dump-map

// works with hamming_search 4.0


#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iterator>
#include <vector>
#include <random>
#include <map>
#include <set>
#include <unordered_set>
#include <string>
#include <cstring>
#include <algorithm>
#include <bitset>
#include <unordered_map>
#include <time.h>
#include <chrono>
#include <new>

#define SL 32
#define SIG_POSSIBLE 268435456



using namespace std;

// prototypes
uint64_t encodekmer(const char s[]);
uint32_t encodekmer_bits(const char *s, vector<int> pos);
bool hd(uint64_t x, uint64_t y, int m);

uint64_t file_read( istream & is, vector <char> & buff );
uint64_t count_lines( const vector <char> & buff, int sz );



int main(int argc, char *argv[]) {

    //typedef unsigned long long int uint64_t;

    using namespace std;

    auto start = chrono::steady_clock::now();
    srand(time(NULL));
    uint64_t kmer_count = 0;

    if (argc <= 7)
    {
        printf("-- Arguments supplied are \n");

        //cout << "buffer\n";
        const int SZ = 1024 * 1024;
        vector <char> buff( SZ );
        ifstream ifs( argv[2]);

        while( int cc = file_read( ifs, buff ) )
        {
            kmer_count += count_lines( buff, cc );
        }

        ifs.close();

    }
    else if (argc > 7)
    {
        printf("Too many arguments supplied.\n");
        exit(0);
    }



    string input_fin = argv[2];
    ifstream fin(input_fin);


    uint32_t p = atoi(argv[4]);
    //int SL = 32;
    uint32_t L = 10;
    float alpha = 0.95;
    uint32_t K = round(log(1 - (pow((1 - alpha), (1 / float(L))))) / log(1 - float(p) / float(SL)));


    cout << "Filename " << input_fin << endl;
    cout << "k-mer count = " << kmer_count << endl;
    cout << "p = " << p << '\n';
    cout << "SL = " << int(SL) << '\n';
    cout << "L = " << L << '\n';
    cout << "alpha = " << alpha << '\n';
    cout << "Using K = " << K << '\n';


    // generate mask
    int l, k, n, sl;

    vector <vector<int> > positions (L, vector<int> (K));

    srand(time(NULL));

    for (l = 0; l < L; l++) {
        vector<int> rand_num;


        for (k = 0; k < K; k++) {
            n = rand() % int(SL);

            if (count(rand_num.begin(), rand_num.end(), n)) {
                k -= 1;
            }
            else {
                rand_num.push_back(n);
            }

        }

        sort(rand_num.begin(), rand_num.end());
        for (int j =0; j < K; j++)
        {
            positions[l][j] = rand_num[j];
            //cout << positions[l][j] << ", ";
        }
        //cout << endl;

        rand_num.clear();
    }


    // write map to file
    string map_name = argv[6];

    FILE *wf;
    wf = fopen(argv[6], "wb");

    if (!wf) {
        cout << "Cannot open file!" << endl;
        return 1;
    }


    fwrite(&p, sizeof(uint32_t), 1, wf);

//    int SL_val = SL;
//    wf.write((char *)&SL_val, sizeof(int8_t));

    fwrite(&L, sizeof(uint32_t), 1, wf);
    fwrite(&alpha, sizeof(float), 1, wf);
    fwrite(&K, sizeof(uint32_t), 1, wf);
    fwrite(&kmer_count, sizeof(uint64_t), 1, wf);
    //cout << "kmer_count " << kmer_count << endl;



    for (l=0; l < L; l++)
    {
        for (int j =0; j < K; j++)
        {
            fwrite(&positions[l][j], sizeof(uint32_t), 1, wf);
            //cout << "l: " << l << " : " << positions[l][j] << endl;
        }

    }




    string line;
    uint64_t li = 0;
    uint64_t b;
    uint32_t sig_hash;
    //uint64_t kmers_processed = 0;

    // write sigs array
    uint64_t total_members_written = 0;
    uint64_t total_sigs_written = 0;


    while (std::getline(fin, line)) {
        if (line.rfind('>', 0) != 0) {
            const char *cline = line.c_str();

            b = encodekmer(cline);
            //cout << b << endl;
            total_members_written += fwrite(&b, sizeof(uint64_t), 1, wf);

            //cout << encode_arr[li] << endl;


            if (li % 1000000 == 0) {
                cout << "-- Encoding  " << li << endl;
            }

            for (l = 0; l < L; l++) {

                sig_hash = encodekmer_bits(cline, positions[l]);
                //cout << sig_hash<< endl;
                total_sigs_written += fwrite(&sig_hash, sizeof(uint32_t), 1, wf);

                //cout <<  *(sigs_arr + l*kmer_count + li) << endl;
            }

            //cout << "Sigs written  " << li << " : " << total_sigs_written << endl;
            //total_sigs_written = 0;


            li += 1;

        }
    }

    cout << "Encodings written " << ": " << total_members_written << endl;

    fin.close();
    fclose(wf);

    auto end = chrono::steady_clock::now();
    cout << "-- Done hashing and writing.  Now writing. Time so far: "
         << chrono::duration_cast<chrono::seconds>(end - start).count() << " seconds" << endl;

    return 0;

}


// function definition
uint64_t encodekmer(const char *s){
    using namespace std;
    uint64_t d = 0;

    for (int i = 0; i <int(SL); i++)
    {
        d = d << 2;

        if (s[i] == 'T'){
            d +=3;
        }

        else if (s[i] == 'G'){
            d +=2;
        }

        else if (s[i] == 'C'){
            d +=1;
        }

        else{
            d +=0;
        }

//        bitset<32> x(d);
//        cout << x << endl;
    }

    return d;
}

bool hd(uint64_t x, uint64_t y, int m){
    uint64_t z = x ^ y;
    int ans = 0;
    //int i = 0;
    for (int i = 0; i < int(SL); i++)
    {
        if (z % 4 != 0)
        {
            ans += 1;
        }

        z = z >> 2;

        if (ans > m)
        {
            return false;
        }
    }

    return true;
}

uint32_t encodekmer_bits(const char *s, vector<int> pos)
{

    using namespace std;
    uint32_t d = 0;

    for (int i = 0; i < pos.size(); i++)
    {
        d = d << 2;

        if (s[pos[i]] == 'T'){
            d +=3;
        }

        else if (s[pos[i]] == 'G'){
            d +=2;
        }

        else if (s[pos[i]] == 'C'){
            d +=1;
        }

        else{
            d +=0;
        }

    }

    return d;
}

uint64_t file_read( istream & is, vector <char> & buff ) {
    is.read( &buff[0], buff.size() );
    return is.gcount();
}

uint64_t count_lines( const vector <char> & buff, int sz ) {
    uint64_t newlines = 0;
    const char * p = &buff[0];
    for ( int i = 0; i < sz; i++ ) {
        if ( p[i] == '>' ) {
            newlines++;
        }
    }
    return newlines;
}

