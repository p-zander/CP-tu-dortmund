#include <fstream>

int main()
{
	auto lcg_init = []( int64_t r0, int64_t a0, int64_t c0, int64_t m0){
		return [=]() mutable{
			r0 = (a0 * r0 + c0) % m0;
			return static_cast<double>(r0) / m0;
		};
	};

	auto lcg1 = lcg_init(1234, 20, 120, 6075);
	auto lcg2 = lcg_init(1234, 137, 187, 256);
	auto lcg3 = lcg_init(123456789, 65539, 0, 2147483648 );
	auto lcg4 = lcg_init(1234, 16807, 0, 2147483647);

	std::ofstream file("lcg.txt");
	file << "#lcg1\tlcg2\tlcg3\tlcg4" << std::endl;
	for (int i=0;i<1e5;i++)
		file << lcg1() << "\t" << lcg2() << "\t" << lcg3() << "\t" << lcg4() <<"\n";
	file.close();
	
	return 0;
}