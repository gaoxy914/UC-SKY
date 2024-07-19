#include "UCSKY.h"

int main(int argc, char const *argv[]) {
    /* code */
    /* int dim = 2;
    int n = 19664;
    double center = 0.5;
    string data_path = "data/iip.dat";
    UCSKY_Solver solver(dim, n, center, 4);
    // solver.gen_ind_data();
    // solver.write_data(data_path.c_str());
    solver.load_data(data_path.c_str());
    // solver.print_data();
    // solver.construct_dg();
    // solver.print_data();
    // solver.loop_bsl();
    // solver.basic_dp();
    // solver.group_dp();
    // solver.branch_bound();
    // solver.sample_fim(1000);
    // solver.greedy();
    // solver.check_prob();
    // solver.construct_dg();
    solver.alpha_greedy(0.002);
    solver.print_ucsky(); */

    if (strcmp(argv[1], "gen-data") == 0) {
        for (int d = 2; d <= 10; d += 2) {
            for (double alpha = 0.2; alpha <= 0.8; alpha += 0.15) {
                UCSKY_Solver solver1(d, 2000000, alpha, 0);
                solver1.gen_ind_data();
                solver1.write_data("data/ind/");
                UCSKY_Solver solver2(d, 2000000, alpha, 0);
                solver2.gen_anti_data();
                solver2.write_data("data/anti/");
                UCSKY_Solver solver3(d, 2000000, alpha, 0);
                solver3.gen_corr_data();
                solver2.write_data("data/corr/");
            }
        }
    } else {
        string data_path = "data/";
        data_path += string(argv[2]);
        int n = atoi(argv[3]);
        int dim = atoi(argv[4]);
        double alpha = atof(argv[5]);
        int l = atoi(argv[6]);
        UCSKY_Solver solver(dim, n, alpha, l);
        solver.load_data(data_path.c_str());

        long long seconds, useconds;
        struct timeval start, end;
        if (strcmp(argv[1], "DP") == 0) {
            gettimeofday(&start, nullptr);
            solver.basic_dp();
            gettimeofday(&end, nullptr);
            seconds = end.tv_sec - start.tv_sec;
            useconds = end.tv_usec - start.tv_usec;
            cout << "-------------------------------" << endl;
            cout << "n = " << n << " d = " << dim << " alpha = " << alpha << " l = " << l << endl;
            cout << argv[1] << " time = " << seconds*1000000 + useconds << endl;
            solver.print_ucsky();
            cout << "-------------------------------" << endl;
        } else if (strcmp(argv[1], "DP+") == 0) {
            gettimeofday(&start, nullptr);
            solver.group_dp();
            gettimeofday(&end, nullptr);
            seconds = end.tv_sec - start.tv_sec;
            useconds = end.tv_usec - start.tv_usec;
            cout << "-------------------------------" << endl;
            cout << "n = " << n << " d = " << dim << " alpha = " << alpha << " l = " << l << endl;
            cout << argv[1] << " time = " << seconds*1000000 + useconds << endl;
            solver.print_ucsky();
            cout << "-------------------------------" << endl;
        } else if (strcmp(argv[1], "GDY") == 0) {
            double step = atof(argv[7]);
            gettimeofday(&start, nullptr);
            solver.alpha_greedy(step);
            gettimeofday(&end, nullptr);
            seconds = end.tv_sec - start.tv_sec;
            useconds = end.tv_usec - start.tv_usec;
            cout << "-------------------------------" << endl;
            cout << "n = " << n << " d = " << dim << " alpha = " << alpha << " l = " << l << " step = " << step << endl;
            cout << argv[1] << " time = " << seconds*1000000 + useconds << endl;
            solver.print_ucsky();
            cout << "-------------------------------" << endl;
        } else if (strcmp(argv[1], "PDBB") == 0) {
            double step = atof(argv[7]);
            gettimeofday(&start, nullptr);
            solver.branch_bound(step);
            gettimeofday(&end, nullptr);
            seconds = end.tv_sec - start.tv_sec;
            useconds = end.tv_usec - start.tv_usec;
            cout << "-------------------------------" << endl;
            cout << "n = " << n << " d = " << dim << " alpha = " << alpha << " l = " << l << " step = " << step << endl;
            cout << argv[1] << " time = " << seconds*1000000 + useconds << endl;
            solver.print_ucsky();
            cout << "-------------------------------" << endl;
        } else if (strcmp(argv[1], "PRE") == 0) {
            int theta = atoi(argv[7]);
            gettimeofday(&start, nullptr);
            solver.sample_fim(theta);
            gettimeofday(&end, nullptr);
            seconds = end.tv_sec - start.tv_sec;
            useconds = end.tv_usec - start.tv_usec;
            cout << "-------------------------------" << endl;
            cout << "n = " << n << " d = " << dim << " alpha = " << alpha << " l = " << l << " theta = " << theta << endl;
            cout << argv[1] << " time = " << seconds*1000000 + useconds << endl;
            solver.print_ucsky();
            cout << "-------------------------------" << endl;
        } else if (strcmp(argv[1], "BSL") == 0) {
            gettimeofday(&start, nullptr);
            solver.loop_bsl();
            gettimeofday(&end, nullptr);
            seconds = end.tv_sec - start.tv_sec;
            useconds = end.tv_usec - start.tv_usec;
            cout << "-------------------------------" << endl;
            cout << "n = " << n << " d = " << dim << " alpha = " << alpha << " l = " << l << endl;
            cout << argv[1] << " time = " << seconds*1000000 + useconds << endl;
            solver.print_ucsky();
            cout << "-------------------------------" << endl;
        } else {
            cout << "wrong alg.\n";
        }
    }
    return 0;
}
