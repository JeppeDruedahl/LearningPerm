// control how much to print in logs
#define LOG_TEST_LEVEL 99
#define LOG_SOLVE_LEVEL 0
#define LOG_SIMULATE_LEVEL 0
#define LOG_CALC_COVS 0
#define LOG_CALC_LOGDIFFS 0

namespace logs {

void test(int print_level, const char * txt, ... )
{

  #if LOG_TEST_LEVEL > 0

    // a. determine behavior
    FILE* log_file;
    if(print_level == 0){ // clean
        log_file = fopen("log_test.txt","w");
    } else if(print_level < LOG_TEST_LEVEL){ // append
        log_file = fopen("log_test.txt","a");
        for(int j = 1; j < print_level; j++){ // indentation
            fprintf(log_file," ");
        }
    } else { // nothing
        return;
    }

    // b. print
    va_list args;
    va_start (args, txt);
    vfprintf (log_file, txt, args);

    // c. close down
    fclose(log_file);
    va_end (args);

  #endif

}

void solve(int print_level, const char * txt, ... )
{

  #if LOG_SOLVE_LEVEL > 0

    // a. determine behavior
    FILE* log_file;
    if(print_level == 0){ // clean
        log_file = fopen("log_solve.txt","w");
    } else if(print_level < LOG_SOLVE_LEVEL){ // append
        log_file = fopen("log_solve.txt","a");
        for(int j = 1; j < print_level; j++){ // indentation
            fprintf(log_file," ");
        }
    } else { // nothing
        return;
    }

    // b. print
    va_list args;
    va_start (args, txt);
    vfprintf (log_file, txt, args);

    // c. close down
    fclose(log_file);
    va_end (args);

  #endif

}

void simulate(int print_level, const char * txt, ... )
{

  #if LOG_SIMULATE_LEVEL > 0

    // a. determine behavior
    FILE* log_file;
    if(print_level == 0){ // clean
        log_file = fopen("log_simulate.txt","w");
    } else if(print_level < LOG_SIMULATE_LEVEL){ // append
        log_file = fopen("log_simulate.txt","a");
        for(int j = 1; j < print_level; j++){ // indentation
            fprintf(log_file," ");
        }
    } else { // nothing
        return;
    }

    // b. print
    va_list args;
    va_start (args, txt);
    vfprintf (log_file, txt, args);

    // c. close down
    fclose(log_file);
    va_end (args);

  #endif

}

void calc_covs(int print_level, const char * txt, ... )
{

  #if LOG_CALC_COVS > 0

    // a. determine behavior
    FILE* log_file;
    if(print_level == 0){ // clean
        log_file = fopen("log_calc_covs.txt","w");
    } else if(print_level < LOG_SIMULATE_LEVEL){ // append
        log_file = fopen("log_calc_covs.txt","a");
        for(int j = 1; j < print_level; j++){ // indentation
            fprintf(log_file," ");
        }
    } else { // nothing
        return;
    }

    // b. print
    va_list args;
    va_start (args, txt);
    vfprintf (log_file, txt, args);

    // c. close down
    fclose(log_file);
    va_end (args);

  #endif

}

void calc_logdiffs(int print_level, const char * txt, ... )
{

  #if LOG_CALC_LOGDIFFS > 0

    // a. determine behavior
    FILE* log_file;
    if(print_level == 0){ // clean
        log_file = fopen("log_calc_logdiffs.txt","w");
    } else if(print_level < LOG_SIMULATE_LEVEL){ // append
        log_file = fopen("log_calc_logdiffs.txt","a");
        for(int j = 1; j < print_level; j++){ // indentation
            fprintf(log_file," ");
        }
    } else { // nothing
        return;
    }

    // b. print
    va_list args;
    va_start (args, txt);
    vfprintf (log_file, txt, args);

    // c. close down
    fclose(log_file);
    va_end (args);

  #endif

}
} // namespace