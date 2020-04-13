// main.cpp -- defining, prototyping, and calling a function

//to compile:
//g++ main_search_v3.4cpp -std = c++11 -o main_search

//to run:
// ./main_search -i G000399765_fna_neighborhood_k32C-map -fq /Users/admin/CLionProjects/hamming_search_1.0/excluded_fna_fq_downSmpl10M
// ./main_search -i G000016385_neighbors_copy2.fa-map -fq /Users/admin/CLionProjects/hamming_search_1.0/excluded_fna_fq_downSmpl10M
// ./main_search -i G000307305_dump-map -fq /Users/admin/CLionProjects/hamming_search_1.0/excluded_fna_fq_downSmpl10M
// works with map from main_map_v4.0

//example: https://theboostcpplibraries.com/boost.multiindex

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iterator>
#include <vector>
#include <map>
#include <string>
#include <cstring>
#include <algorithm>
#include <bitset>
#include <unordered_map>
#include <time.h>
#include <dirent.h>
#include <sys/types.h>
#include <chrono>
#include <new>
#include <algorithm>


#include <boost/multi_index_container.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/member.hpp>

#include <boost/multi_index/sequenced_index.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/random_access_index.hpp>



using namespace boost::multi_index;

struct kmer
{
    uint64_t encoding;
    uint32_t sig_l0;
    uint32_t sig_l1;
    uint32_t sig_l2;
    uint32_t sig_l3;
    uint32_t sig_l4;
    uint32_t sig_l5;
    uint32_t sig_l6;
    uint32_t sig_l7;
    uint32_t sig_l8;
    uint32_t sig_l9;

};

typedef multi_index_container<
        kmer,
        indexed_by<
                ordered_unique<
                        member<
                                kmer, uint64_t, &kmer::encoding
                        >
                >,
                ordered_non_unique<
                        member<
                                kmer, uint32_t, &kmer::sig_l0
                        >
                >,
                ordered_non_unique<
                        member<
                                kmer, uint32_t, &kmer::sig_l1
                        >
                >,
                ordered_non_unique<
                        member<
                                kmer, uint32_t, &kmer::sig_l2
                        >
                >,
                ordered_non_unique<
                        member<
                                kmer, uint32_t, &kmer::sig_l3
                        >
                >,
                ordered_non_unique<
                        member<
                                kmer, uint32_t, &kmer::sig_l4
                        >
                >,
                ordered_non_unique<
                        member<
                                kmer, uint32_t, &kmer::sig_l5
                        >
                >,
                ordered_non_unique<
                        member<
                                kmer, uint32_t, &kmer::sig_l6
                        >
                >,
                ordered_non_unique<
                        member<
                                kmer, uint32_t, &kmer::sig_l7
                        >
                >,
                ordered_non_unique<
                        member<
                                kmer, uint32_t, &kmer::sig_l8
                        >
                >,
                ordered_non_unique<
                        member<
                                kmer, uint32_t, &kmer::sig_l9
                        >
                >
        >
> kmer_multi;



#define SL 32




// prototypes
uint64_t encodekmer(const char s[]);
uint64_t encodekmer_rev(const char s[]);
uint32_t encodekmer_bits(const char *s, std::vector<int> pos);
uint32_t encodekmer_bits_rev(const char *s, std::vector<int> pos);
bool hd(uint64_t x, uint64_t y, int m);
std::vector<std::string> list_dir(const char *path);



int main(int argc, char *argv[]) {


    std::cout << "after boost install !!" << std::endl;

    kmer_multi kmers;





    std::cout << "done testing boost" << std::endl;




    using namespace std;

    //typedef unsigned long long int  uint64_t;
    auto start = chrono::steady_clock::now();


    if (argc <= 5) {
        printf("-- Arguments supplied are \n");
        //printf("Filename %s\n", argv[2]);
        //printf("p = %d\n", atoi(argv[4]));
    } else if (argc > 5) {
        printf("Too many arguments supplied.\n");
        exit(0);
    }

    //read map

    char *input_fin = argv[2];

    FILE *f = fopen(input_fin, "rb");

    // read parameters from input file
    uint32_t p;
    fread(&p, sizeof(uint32_t), 1, f);

    uint32_t L;
    fread(&L, sizeof(uint32_t), 1, f);

    float alpha;
    fread(&alpha, sizeof(float), 1, f);

    uint32_t K;
    fread(&K, sizeof(uint32_t), 1, f);

    uint64_t kmer_count;
    fread(&kmer_count, sizeof(uint64_t), 1, f);

    string line;
    int l;


    cout << "Map " << input_fin << endl;
    cout << "k-mer count = " << kmer_count << endl;
    cout << "p = " << p << '\n';
    cout << "SL = " << SL << '\n';
    cout << "L = " << L << '\n';
    cout << "alpha = " << alpha << '\n';
    cout << "Using K = " << K << '\n';



    // read mask
    vector<vector<int> > positions(L, vector<int>(K));

    //vector<vector<u_int32_t > > sigs(L, vector<uint32_t>(kmer_count));
    vector<uint64_t> encods(kmer_count);



    for (l = 0; l < L; l++) {

        for (int j = 0; j < K; j++) {
            fread(&positions[l][j], sizeof(uint32_t), 1, f);
            //cout << "l: " << l << " : " << positions[l][j] << endl;
        }

    }


    uint64_t num_pairs = 0;
    uint64_t sigs_pairs = 0;


    for (int enc = 0; enc < kmer_count; enc++) {
        num_pairs += fread(&encods[enc], sizeof(uint64_t), 1, f);
        //cout << encods[enc] << endl;

        uint32_t s0, s1, s2, s3, s4, s5, s6, s7, s8, s9;


            sigs_pairs+= fread(&s0, sizeof(uint32_t), 1, f);
            sigs_pairs+= fread(&s1, sizeof(uint32_t), 1, f);
            sigs_pairs+= fread(&s2, sizeof(uint32_t), 1, f);
            sigs_pairs+= fread(&s3, sizeof(uint32_t), 1, f);
            sigs_pairs+= fread(&s4, sizeof(uint32_t), 1, f);
            sigs_pairs+= fread(&s5, sizeof(uint32_t), 1, f);
            sigs_pairs+= fread(&s6, sizeof(uint32_t), 1, f);
            sigs_pairs+= fread(&s7, sizeof(uint32_t), 1, f);
            sigs_pairs+= fread(&s8, sizeof(uint32_t), 1, f);
            sigs_pairs+= fread(&s9, sizeof(uint32_t), 1, f);



//        cout << "Sigs read  " << enc << " : " << sigs_pairs << endl;
//        cout << s0 << endl;
//        cout << s1 << endl;
//        cout << s2 << endl;
//        cout << s3 << endl;
//        cout << s4 << endl;
//        cout << s5 << endl;
//        cout << s6 << endl;
//        cout << s7 << endl;
//        cout << s8 << endl;
//        cout << s9 << endl;

        kmers.insert({encods[enc], s0, s1, s2, s3, s4, s5, s6, s7, s8, s9});
        //sigs_pairs = 0;
    }

    cout << "Sigs read  " << " : " << sigs_pairs/L << endl;
    cout << "Encodings read " << " : " << num_pairs << endl;



    fclose(f);


    auto end = chrono::steady_clock::now();
    cout << "-- Done reading. Now matching. Time so far: "
         << chrono::duration_cast<chrono::seconds>(end - start).count() << " seconds" << endl;



    // read input fastq
    const char *dir = argv[4];

    vector<string> file_list;
    if (dir == NULL) {
        return (1);
    } else {
        file_list = list_dir(dir);
    }

    int file_count = file_list.size();
    for (int f = 0; f < file_count; f++) {

        string input_fq = file_list[f];
        cout << input_fq << endl;

        stringstream stream;
        stream << std::fixed << std::setprecision(2) << alpha;
        std::string alpha_s = stream.str();

        string input_fq_truct = input_fq.substr(input_fq.find_last_of("/") + 1);

        //string input_map = argv[2];
        //input_map = input_map.substr(0, input_map.find("_"));

        //    string output_fname = "test_k" + to_string(int(SL)) + "_l" + to_string(int(L)) + "_p" +
        //            to_string(int(p))+ "_alpha" + alpha_s + "_" + input_map + "_" +input_fq_truct;

        string output_fname = "test_k" + to_string(int(SL)) + "_l" + to_string(int(L)) + "_p" +
                              to_string(int(p)) + "_alpha" + alpha_s + "_" + input_fq_truct;

        ifstream ifs(input_fq);

        ofstream outputFile;
        outputFile.open(output_fname);

        uint64_t lines_read = 0;
        uint64_t reads_matched = 0;

        string line_of_file;
        string name;

        uint64_t b;
        uint32_t kmer_sig;


        while (!ifs.eof()) {
            getline(ifs, line_of_file);

            if (ifs) {

                if (lines_read % 4 == 0) {
                    name = line_of_file;
                }

                if (lines_read % 4 == 1) {
                    int matched = 0;
                    //cout << lines_read << endl;


                    for (int i = 0; i < line_of_file.length() - int(SL) + 1; i++) {
                        string kmer_str = line_of_file.substr(i, int(SL));
                        //cout << kmer_str << endl;
                        const char *ckmer = kmer_str.c_str();
                        b = encodekmer(ckmer);
                        //cout << b << endl;
                        bool kmerfound = false;


                        for (int funci = 0; funci < L; funci++) {
                            kmer_sig = encodekmer_bits(ckmer, positions[funci]);
                            //cout << kmer_sig << endl;


                            if (funci == 0){

                                //cout << funci << endl;

                                auto &sigl0_index = kmers.get<1>();
                                auto it = sigl0_index.lower_bound(kmer_sig);
                                auto kmer_end = sigl0_index.upper_bound(kmer_sig);
                                for (; it != kmer_end; ++it){
                                    if (hd(b, it->encoding, p)) {
                                        //cout << funci << endl;
                                        //cout << it->encoding << endl;
                                        matched += 1;
                                        kmerfound = true;
                                        break;
                                    }
                                }
                            }

                            else if (funci == 1){
                                auto &sigl1_index = kmers.get<2>();
                                auto it = sigl1_index.lower_bound(kmer_sig);
                                auto kmer_end = sigl1_index.upper_bound(kmer_sig);
                                for (; it != kmer_end; ++it)
                                {
                                    if (hd(b, it->encoding, p)) {
                                        //cout << funci << endl;
                                        matched += 1;
                                        kmerfound = true;
                                        break;
                                    }
                                }

                            }
                            else if (funci == 2){
                                auto &sigl2_index = kmers.get<3>();
                                auto it = sigl2_index.lower_bound(kmer_sig);
                                auto kmer_end = sigl2_index.upper_bound(kmer_sig);
                                for (; it != kmer_end; ++it)
                                {
                                    if (hd(b, it->encoding, p)) {
                                        //cout << funci << endl;
                                        matched += 1;
                                        kmerfound = true;
                                        break;
                                    }
                                }

                            }
                            else if (funci == 3){

                                auto &sigl3_index = kmers.get<4>();
                                auto it = sigl3_index.lower_bound(kmer_sig);
                                auto kmer_end = sigl3_index.upper_bound(kmer_sig);
                                for (; it != kmer_end; ++it)
                                {
                                    if (hd(b, it->encoding, p)) {
                                        //cout << funci << endl;
                                        matched += 1;
                                        kmerfound = true;
                                        break;
                                    }
                                }

                            }
                            else if (funci == 4){
                                auto &sigl4_index = kmers.get<5>();
                                auto it = sigl4_index.lower_bound(kmer_sig);
                                auto kmer_end = sigl4_index.upper_bound(kmer_sig);
                                for (; it != kmer_end; ++it)
                                {
                                    if (hd(b, it->encoding, p)) {
                                        //cout << funci << endl;
                                        matched += 1;
                                        kmerfound = true;
                                        break;
                                    }
                                }


                            }
                            else if (funci == 5){
                                auto &sigl5_index = kmers.get<6>();
                                auto it = sigl5_index.lower_bound(kmer_sig);
                                auto kmer_end = sigl5_index.upper_bound(kmer_sig);
                                for (; it != kmer_end; ++it)
                                {
                                    if (hd(b, it->encoding, p)) {
                                        //cout << funci << endl;
                                        matched += 1;
                                        kmerfound = true;
                                        break;
                                    }
                                }

                            }
                            else if (funci == 6){
                                auto &sigl6_index = kmers.get<7>();
                                auto it = sigl6_index.lower_bound(kmer_sig);
                                auto kmer_end = sigl6_index.upper_bound(kmer_sig);
                                for (; it != kmer_end; ++it)
                                {
                                    if (hd(b, it->encoding, p)) {
                                        //cout << funci << endl;
                                        matched += 1;
                                        kmerfound = true;
                                        break;
                                    }
                                }

                            }
                            else if (funci == 7){
                                auto &sigl7_index = kmers.get<8>();
                                auto it = sigl7_index.lower_bound(kmer_sig);
                                auto kmer_end = sigl7_index.upper_bound(kmer_sig);
                                for (; it != kmer_end; ++it)
                                {
                                    if (hd(b, it->encoding, p)) {
                                        //cout << funci << endl;
                                        matched += 1;
                                        kmerfound = true;
                                        break;
                                    }
                                }

                            }
                            else if (funci == 8){
                                auto &sigl8_index = kmers.get<9>();
                                auto it = sigl8_index.lower_bound(kmer_sig);
                                auto kmer_end = sigl8_index.upper_bound(kmer_sig);
                                for (; it != kmer_end; ++it)
                                {
                                    if (hd(b, it->encoding, p)) {
                                        //cout << funci << endl;
                                        matched += 1;
                                        kmerfound = true;
                                        break;
                                    }
                                }

                            }
                            else if (funci == 9){
                                auto &sigl9_index = kmers.get<10>();
                                auto it = sigl9_index.lower_bound(kmer_sig);
                                auto kmer_end = sigl9_index.upper_bound(kmer_sig);
                                for (; it != kmer_end; ++it)
                                {
                                    if (hd(b, it->encoding, p)) {
                                        //cout << funci << endl;
                                        matched += 1;
                                        kmerfound = true;
                                        break;
                                    }
                                }

                            }

                                if (kmerfound) {
                                    break;
                                }
                            }
                        }



                        // try reverse complement
                        if (matched == 0) {

                            int len = strlen(line_of_file.c_str());
                            char swap;

                            for (int i = 0; i < len / 2; i++) {
                                swap = line_of_file[i];
                                line_of_file[i] = line_of_file[len - i - 1];
                                line_of_file[len - i - 1] = swap;
                            }


                            for (int i = 0; i < line_of_file.length() - int(SL) + 1; i++) {
                                string kmer_str = line_of_file.substr(i, int(SL));
                                const char *ckmer = kmer_str.c_str();
                                b = encodekmer_rev(ckmer);
                                //cout << b << endl;


                                bool kmerfound = false;

                                for (int funci = 0; funci < L; funci++) {

                                    kmer_sig = encodekmer_bits_rev(ckmer, positions[funci]);
                                    //cout << kmer_sig << endl;


                                    if (funci == 0) {

                                        auto &sigl0_index = kmers.get<1>();
                                        auto it = sigl0_index.lower_bound(kmer_sig);
                                        auto kmer_end = sigl0_index.upper_bound(kmer_sig);
                                        for (; it != kmer_end; ++it)
                                        {
                                            if (hd(b, it->encoding, p)) {
                                                //cout << "rev" << funci << endl;
                                                matched += 1;
                                                kmerfound = true;
                                                break;
                                            }
                                        }
                                    } else if (funci == 1){
                                        auto &sigl1_index = kmers.get<2>();
                                        auto it = sigl1_index.lower_bound(kmer_sig);
                                        auto kmer_end = sigl1_index.upper_bound(kmer_sig);
                                        for (; it != kmer_end; ++it)
                                        {
                                            if (hd(b, it->encoding, p)) {
                                                //cout << "rev" << funci << endl;
                                                matched += 1;
                                                kmerfound = true;
                                                break;
                                            }
                                        }

                                    }
                                    else if (funci == 2){
                                        auto &sigl2_index = kmers.get<3>();
                                        auto it = sigl2_index.lower_bound(kmer_sig);
                                        auto kmer_end = sigl2_index.upper_bound(kmer_sig);
                                        for (; it != kmer_end; ++it)
                                        {
                                            if (hd(b, it->encoding, p)) {
                                                //cout << "rev" << funci << endl;
                                                matched += 1;
                                                kmerfound = true;
                                                break;
                                            }
                                        }

                                    }
                                    else if (funci == 3){

                                        auto &sigl3_index = kmers.get<4>();
                                        auto it = sigl3_index.lower_bound(kmer_sig);
                                        auto kmer_end = sigl3_index.upper_bound(kmer_sig);
                                        for (; it != kmer_end; ++it)
                                        {
                                            if (hd(b, it->encoding, p)) {
                                                //cout << "rev" << funci << endl;
                                                matched += 1;
                                                kmerfound = true;
                                                break;
                                            }
                                        }

                                    }
                                    else if (funci == 4){
                                        auto &sigl4_index = kmers.get<5>();
                                        auto it = sigl4_index.lower_bound(kmer_sig);
                                        auto kmer_end = sigl4_index.upper_bound(kmer_sig);
                                        for (; it != kmer_end; ++it)
                                        {
                                            if (hd(b, it->encoding, p)) {
                                                //cout << "rev" << funci << endl;
                                                matched += 1;
                                                kmerfound = true;
                                                break;
                                            }
                                        }


                                    }
                                    else if (funci == 5){
                                        auto &sigl5_index = kmers.get<6>();
                                        auto it = sigl5_index.lower_bound(kmer_sig);
                                        auto kmer_end = sigl5_index.upper_bound(kmer_sig);
                                        for (; it != kmer_end; ++it)
                                        {
                                            if (hd(b, it->encoding, p)) {
                                                //cout << "rev" << funci << endl;
                                                matched += 1;
                                                kmerfound = true;
                                                break;
                                            }
                                        }

                                    }
                                    else if (funci == 6){
                                        auto &sigl6_index = kmers.get<7>();
                                        auto it = sigl6_index.lower_bound(kmer_sig);
                                        auto kmer_end = sigl6_index.upper_bound(kmer_sig);
                                        for (; it != kmer_end; ++it)
                                        {
                                            if (hd(b, it->encoding, p)) {
                                                //cout << "rev" << funci << endl;
                                                matched += 1;
                                                kmerfound = true;
                                                break;
                                            }
                                        }

                                    }
                                    else if (funci == 7){
                                        auto &sigl7_index = kmers.get<8>();
                                        auto it = sigl7_index.lower_bound(kmer_sig);
                                        auto kmer_end = sigl7_index.upper_bound(kmer_sig);
                                        for (; it != kmer_end; ++it)
                                        {
                                            if (hd(b, it->encoding, p)) {
                                                //cout << "rev" << funci << endl;
                                                matched += 1;
                                                kmerfound = true;
                                                break;
                                            }
                                        }

                                    }
                                    else if (funci == 8){
                                        auto &sigl8_index = kmers.get<9>();
                                        auto it = sigl8_index.lower_bound(kmer_sig);
                                        auto kmer_end = sigl8_index.upper_bound(kmer_sig);
                                        for (; it != kmer_end; ++it)
                                        {
                                            if (hd(b, it->encoding, p)) {
                                                //cout << "rev" << funci << endl;
                                                matched += 1;
                                                kmerfound = true;
                                                break;
                                            }
                                        }

                                    }
                                    else if (funci == 9){
                                        auto &sigl9_index = kmers.get<10>();
                                        auto it = sigl9_index.lower_bound(kmer_sig);
                                        auto kmer_end = sigl9_index.upper_bound(kmer_sig);
                                        for (; it != kmer_end; ++it)
                                        {
                                            if (hd(b, it->encoding, p)) {
                                                //cout << "rev" << funci << endl;
                                                matched += 1;
                                                kmerfound = true;
                                                break;
                                            }
                                        }

                                    }

                                    if (kmerfound) {
                                        break;
                                    }
                                }




                            }

                        }



                        if (matched > 0) {
                            outputFile << ">" << name << " " << matched << endl;
                            outputFile << line_of_file << endl;
                            reads_matched += 1;
                        }
//                        else{
//                        cout << lines_read << " : " << line_of_file;
//                    }

                    }

                    if (lines_read % 100000 == 0) {
                        cout << lines_read << " " << reads_matched << endl;
                    }

                    ++lines_read;
                } else break;
            }

            cout << lines_read << " " << reads_matched << endl;
            end = chrono::steady_clock::now();
            cout << "-- Done matching. Time so far: " << chrono::duration_cast<chrono::seconds>(end - start).count()
                 << " seconds" << endl;

            ifs.close();
            outputFile.close();
        }


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

    }

    return d;
}


uint64_t encodekmer_rev(const char *s){
    using namespace std;
    uint64_t d = 0;

    for (int i = 0; i <int(SL); i++)
    {
        d = d << 2;

        if (s[i] == 'A'){
            d +=3;
        }

        else if (s[i] == 'C'){
            d +=2;
        }

        else if (s[i] == 'G'){
            d +=1;
        }

        else{
            d +=0;
        }

    }

    return d;
}




bool hd(uint64_t x, uint64_t y, int m){
    using namespace std;

    uint64_t z = x ^ y;
    int ans = 0;

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


std::vector<std::string> list_dir(const char *path) {
    using namespace std;
    vector<string> userString;
    struct dirent *entry;
    DIR *dir = opendir(path);

    while ((entry = readdir(dir)) != NULL) {
        string temp = string(path) + "/" + entry->d_name;
        userString.push_back(temp);
    }
    closedir(dir);

    return(userString);
}

uint32_t encodekmer_bits(const char *s, std::vector<int> pos)
{

    using namespace std;
    uint32_t d = 0;

    for (int i = 0; i < pos.size(); i++)
    {
        d = d << 2;
        //cout << pos.size() << endl;

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

uint32_t encodekmer_bits_rev(const char *s, std::vector<int> pos){
    using namespace std;
    uint32_t d = 0;

    for (int i = 0; i < pos.size(); i++)
    {
        d = d << 2;

        if (s[pos[i]] == 'A'){
            d +=3;
        }

        else if (s[pos[i]] == 'C'){
            d +=2;
        }

        else if (s[pos[i]] == 'G'){
            d +=1;
        }

        else{
            d +=0;
        }

    }

    return d;
}