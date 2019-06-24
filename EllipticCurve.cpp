#include "EllipticCurve.hpp"
#include <cstdint>
#include <iostream>
#include <chrono>
#include <string>

/// Driver program for the elliptic curve impl
int main()
{
	// Benchmark timestamp
	auto now_ts = std::chrono::system_clock::now();

	// Parametrized elliptic curve
	EllipticCurve<uint64_t> curve(10374515419998325367ull, 2, 4);
	
	// Get a point that exists regardless of GF
	auto point = curve(2, 4);

	// Print
	std::cout << point << '\n';

	// Perform multiplication
	point *= uint64_t(16989515419998325367ull);

	// Print the multiplied point
	std::cout << point << '\n';

	// Print the inverse of the point
	auto inverse_point = point.Inverse();
	std::cout << inverse_point << '\n';

	// Print the order of the first point
	std::cout << point.Order() << '\n';

	// Print the order of the inverse of the first point
	std::cout << inverse_point.Order() << '\n';

	// Benchmark finished - get the timestamp
	auto end_ts = std::chrono::system_clock::now();

	// Output the difference
	std::cout << std::chrono::duration_cast<std::chrono::microseconds>(end_ts - now_ts).count() << '\n';
	return 0;
}