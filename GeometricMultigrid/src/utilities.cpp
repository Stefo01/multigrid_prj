#include "allIncludes.hpp"

void Utils::Initialization_for_N(int argc, char** argv, size_t &N, double &alpha, double &width, int &level, int &functions_to_test, SMOOTHERS &sm){
    std::string nn = "-n";
    std::string aa = "-a";
    std::string ww = "-w";
    std::string lv = "-ml";
    std::string smooth = "-smt";
    std::string help = "--help";
    std::string fti = "-test";
    alpha = DEFAULT_ALPHA;
    N = DEFAULT_N;
    width = DEFAULT_WIDTH;
    level = DEFAULT_LEVEL;
    functions_to_test = DEFAULT_TEST;
    sm = DEFAULT_METHOD;

    if(argc < 2)
    {
        std::cout<<"Inserted by default N = "<<DEFAULT_N<<std::endl;
        std::cout<<"Inserted by default alpha = "<<DEFAULT_ALPHA<<std::endl;
        std::cout<<"Inserted by default width = "<<DEFAULT_WIDTH<<std::endl;
        std::cout<<"Inserted by default multigrid level = "<<DEFAULT_LEVEL<<std::endl;
        std::cout<<"Inserted by default test number "<<DEFAULT_TEST<<std::endl;
        std::cout<<"Inserted by default Smooter number "<<DEFAULT_METHOD<<std::endl;
    } 
    else{
        for(int i= 0; i < argc; i++)
        {
            if(nn.compare(argv[i]) == 0 && ( i+1 < argc ) )
            {
                try{
                    N = std::stoi(argv[i+1]);
                    std::cout<<"Inserted N = "<<N<<std::endl;
                    if(N <= 0)
                    {
                        std::cout<<"Error: Please, insert a valid N value"<<std::endl;
                        std::exit(1);
                    }
                }catch(std::exception&){
                    std::cout<<"Error: Please, insert a number after -n"<<std::endl;
                    std::exit(1);
                }
                
            }
            else if(aa.compare(argv[i]) == 0 && ( i+1 < argc ) )
            {
                try{
                    alpha = std::atof(argv[i+1]);
                    std::cout<<"Inserted alpha = "<<alpha<<std::endl;
                }catch(std::exception&){
                    std::cout<<"Error: Please, insert a double after -a"<<std::endl;
                    std::exit(1);
                }
            }
            else if(lv.compare(argv[i]) == 0 && ( i+1 < argc ) )
            {
                try{
                    level = std::stoi(argv[i+1]);
                    std::cout<<"Inserted level = "<<level<<std::endl;
                    if(level <= 0)
                    {
                        std::cout<<"Error: Please, insert a valid level"<<std::endl;
                        std::exit(1);
                    }
                }catch(std::exception&){
                    std::cout<<"Error: Please, insert a number after -ml"<<std::endl;
                    std::exit(1);
                }
            }
            else if(smooth.compare(argv[i]) == 0 && ( i+1 < argc ) )
            {
                try{
                    sm = static_cast<SMOOTHERS>(std::stoi(argv[i+1]));
                    std::cout<<"Inserted Smoother number = "<<sm<<std::endl;
                    if(sm >= SMOOTHERS::SMOOTHERS_END)
                        sm = DEFAULT_METHOD;
                }catch(std::exception&){
                    std::cout<<"Error: Please, insert a number after -smt"<<std::endl;
                    std::exit(1);
                }
            }
            else if(fti.compare(argv[i]) == 0 && ( i+1 < argc ) )
            {
                try{
                    functions_to_test = std::stoi(argv[i+1]);
                    std::cout<<"Inserted test number = "<<functions_to_test<<std::endl;
                    if(functions_to_test < 0)
                    {
                        std::cout<<"Error: Please, insert a valid test number"<<std::endl;
                        std::exit(1);
                    }
                }catch(std::exception&){
                    std::cout<<"Error: Please, insert a double after -test"<<std::endl;
                    std::exit(1);
                }
            }
            else if(ww.compare(argv[i]) == 0 && ( i+1 < argc ) )
            {
                try{
                    width = std::atof(argv[i+1]);
                    std::cout<<"Inserted width = "<<width<<std::endl;
                    if(width <= 0)
                    {
                        std::cout<<"Error: Please, insert a valid width"<<std::endl;
                        std::exit(1);
                    }
                }catch(std::exception&){
                    std::cout<<"Error: Please, insert a double after -w"<<std::endl;
                    std::exit(1);
                }
            }
            else if(help.compare(argv[i]) == 0)
            {
                std::cout << "Usage: ./Multigrid [OPTIONS]\n"<<std::endl
                << "Options:"<<std::endl
                << "  -n, insert number of spaces"<<std::endl
                << "  -a, specifies differential constant"<<std::endl
                << "  -w, insert the Width of the rectangle domain"<<std::endl
                << "  -ml, insert multigrid level"<<std::endl
                << "  -test, insert type of function in input to test it"<<std::endl
                << "  -smt, you can choose your favourite smoother"<<std::endl
                << "  --help, Display this help message"<<std::endl;
                std::exit(1);
            }
            else if(nn.compare(argv[i]) == 0 && ( i+1 == argc )) {
                std::cout<<"Error: Please, insert something"<<std::endl;
                exit(1);
            }
        } 
    } 
}


void Utils::init_test_functions(std::function<double(const double, const double)> &f, std::function<double(const double, const double)> &g, int i)
{
    using function_to_c = std::function<double(const double, const double)>;
    std::array< function_to_c,  6> functions_to_choose_to_test ={
        // FFF
        [] (const double x, const double y) { return 1.; },
        [] (const double x, const double y) { return 0.; },
        [] (const double x, const double y) { return -5.0 * exp(x) * exp(-2.0 * y); }, 
        [] (const double x, const double y) { return exp(x) * exp(-2.0 * y); }, 
        [] (const double x, const double y) { double r = std::sqrt(x * x + y * y); return r != 0.0 ? -30. * (std::cos(30. * r) / r - 30. * std::sin(30. * r)) : 0.0; }, 
        [] (const double x, const double y) { return std::sin(30. * std::sqrt(x * x + y * y)); }
        // END
    };  

    size_t max_index = functions_to_choose_to_test.size() / 2 - 1;
    if (i < 0 || static_cast<size_t>(i) > max_index) {
        f = functions_to_choose_to_test[0];
        g = functions_to_choose_to_test[1];
        std::cout << "Warning: Invalid test case index. Default test case selected.\n";
    } else {
        f = functions_to_choose_to_test[i * 2];
        g = functions_to_choose_to_test[i * 2 + 1];
    }

}
