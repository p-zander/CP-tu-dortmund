#include <assert.h>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <tuple>

int main()
{
	auto lcg_init = []( int64_t r0, int64_t a0, int64_t c0, int64_t m0){
		return [=]() mutable{
			r0 = (a0 * r0 + c0) % m0;
			return static_cast<double>(r0) / m0;
		};
	};

	auto lcg1 = lcg_init(123456789, 65539, 0, 2147483648 );
	auto lcg2 = lcg_init(123456789, 65539, 0, 2147483648 );
	auto lcg3 = lcg_init(123456789, 65539, 0, 2147483648 );
	auto lcg4 = lcg_init(123456789, 65539, 0, 2147483648 );

	auto bm = [&](){
		auto x1 = lcg1();
		auto x2 = lcg1();
		auto y1 = sqrt(-2 * log(x1)) * cos(2 * M_PI * x2);
		auto y2 = sqrt(-2 * log(x1)) * sin(2 * M_PI * x2);
		return std::make_tuple(y1,y2);
	};

	std::function<double(unsigned int)> zgs;
    zgs = [&](unsigned int n){
	    assert(n!=0);
        if (n == 1){
            return (lcg2() - 0.5) / sqrt(static_cast<double>(n)/1);
        }else
            return (lcg2() - 0.5) / sqrt(static_cast<double>(n)/1) + zgs(n-1);
    };

    auto inv = [&](){return pow(lcg3(),1.0/3);};

	std::function<double()> rueck;
    rueck = [&](){
		auto x1 = lcg4();
		auto x2 = lcg4();
		if (x2*M_PI<(sin(x1*M_PI)/2)) return M_PI * x1;
		else return rueck();
    };



	double x2,x1;
	std::ofstream file("rand.txt");
	// file << "#bm\tlcg2\tlcg3\tlcg4" << std::endl;
	for (int i=0;i<1e5;i++){
		std::tie(x1,x2) = bm();
		file << x1 << "\t" << zgs(5000) << "\t" << rueck() << "\t" << inv() << std::endl;
		file << x2 << "\t" << zgs(5000) << "\t" << rueck() << "\t" << inv() << std::endl;
	}
	file.close();


	return 0;
}