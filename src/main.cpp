#include "UCSKY.h"

int main(int argc, char const *argv[]) {
    /* code */
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
                solver3.write_data("data/corr/");
            }
        }
    } else if (strcmp(argv[1], "gen-data-2d") == 0) {
        int m = atoi(argv[2]);
        UCSKY_Solver solver(2, 100000, 0.5, 0);
        solver.gen_ind_data_2d(m);
        solver.write_data("data/ind/");
    } else {
        string data_path = "data/";
        data_path += string(argv[2]);
        int n = atoi(argv[3]);
        int dim = atoi(argv[4]);
        double alpha = atof(argv[5]);
        int l = atoi(argv[6]);
        UCSKY_Solver solver(dim, n, alpha, l);
        solver.load_data(data_path.c_str());
        double step = 0.01;
        int theta = 160;
        long long seconds, useconds;
        struct timeval start, end;
        int nR = 0;
        if (strcmp(argv[1], "DP") == 0) {
            gettimeofday(&start, nullptr);
            // nR = solver.reduce_data(step);
            // solver.basic_dp();
            solver.basic_dp2();
            // solver.basic_dp3();
            gettimeofday(&end, nullptr);
            seconds = end.tv_sec - start.tv_sec;
            useconds = end.tv_usec - start.tv_usec;
        } else if (strcmp(argv[1], "DP+") == 0) {
            gettimeofday(&start, nullptr);
            // nR = solver.reduce_data(step);
            // solver.group_dp();
            // solver.group_dp2();
            solver.group_dp3();
            gettimeofday(&end, nullptr);
            seconds = end.tv_sec - start.tv_sec;
            useconds = end.tv_usec - start.tv_usec;
        } else if (strcmp(argv[1], "GDY") == 0) {
            step = atof(argv[7]);
            gettimeofday(&start, nullptr);
            solver.greedy_plus(step);
            gettimeofday(&end, nullptr);
            seconds = end.tv_sec - start.tv_sec;
            useconds = end.tv_usec - start.tv_usec;
        } else if (strcmp(argv[1], "PDBB") == 0) {
            gettimeofday(&start, nullptr);
            nR = solver.reduce_data(step);
            // cout << nR << endl;
            solver.branch_bound();
            gettimeofday(&end, nullptr);
            seconds = end.tv_sec - start.tv_sec;
            useconds = end.tv_usec - start.tv_usec;
        } else if (strcmp(argv[1], "PRE") == 0) {
            theta = atoi(argv[7]);
            gettimeofday(&start, nullptr);
            solver.sample_FIM2(theta);
            gettimeofday(&end, nullptr);
            seconds = end.tv_sec - start.tv_sec;
            useconds = end.tv_usec - start.tv_usec;
        } else if (strcmp(argv[1], "SKY") == 0) {
            gettimeofday(&start, nullptr);
            solver.sky_2n();
            gettimeofday(&end, nullptr);
            seconds = end.tv_sec - start.tv_sec;
            useconds = end.tv_usec - start.tv_usec;
        } else if (strcmp(argv[1], "BSL") == 0) {
            gettimeofday(&start, nullptr);
            nR = solver.reduce_data_bsl();
            // cout << nR << endl;
            solver.loop_bsl();
            gettimeofday(&end, nullptr);
            seconds = end.tv_sec - start.tv_sec;
            useconds = end.tv_usec - start.tv_usec;
        } else if (strcmp(argv[1], "test") == 0) {
            // nR = solver.reduce_data(step);
        } else {
            cout << "wrong alg.\n";
        }
        cout << "-------------------------------" << endl;
        cout << data_path << " n = " << n << " d = " << dim << " alpha = " << alpha << " l = " << l << " step = " << step << " theta = " << theta << endl;
        cout << argv[1] << " time = " << seconds*1000000 + useconds << " nR = " << nR << endl;
        // solver.print_ucsky();
        solver.check_ucsky();
        cout << "-------------------------------" << endl;
    }
    return 0;
}
