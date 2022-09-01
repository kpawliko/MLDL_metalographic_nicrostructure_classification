/* benchmark.cpp
 *
 * A complex tool that exposes most functionalities and
 * allows for benchmarking methods in different conditions.
 *
 * Parameters are passed in commandline in form:
 * name=value
 *
 * Running the program without any arguments will print available options.
 */
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <numeric>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <filesystem>

#include "timer.hpp"
#include "papi.hpp"
#include "random.hpp"
#include "rca.hpp"
#include "visualize.hpp"
#include "argument_parser.hpp"
#include "get_simulation.hpp"

using namespace std;
using namespace rca;

using S = rca::Grain;

extern template class rca::Material<float, S, PeriodicDisabled>;
extern template class rca::Material<double, S, PeriodicDisabled>;
extern template class rca::Material<long double, S, PeriodicDisabled>;

extern template class rca::Material<float, S, PeriodicEnabled>;
extern template class rca::Material<double, S, PeriodicEnabled>;
extern template class rca::Material<long double, S, PeriodicEnabled>;

extern template class rca::MaterialMPI<float, rca::Grain>;
extern template class rca::MaterialMPI<double, rca::Grain>;
extern template class rca::MaterialMPI<long double, rca::Grain>;

std::vector<std::string> split(std::string s, std::string delim) {
	std::vector<std::string> fields;
	size_t pos;
	while ((pos = s.find(delim)) != std::string::npos) {
		auto field = s.substr(0, pos);
		fields.push_back(field);
		s.erase(0, pos + delim.length());
	}
	fields.push_back(s);
	return fields;
}

void setOmpSchedule(ArgumentParser &ap) {
	auto schedule_name = ap.get<std::string>("omp-schedule", "static");
	auto chunk_size = ap.get<int>("omp-chunk-size", 1);

	omp_sched_t sched;
	if (schedule_name == "static") {
		sched = omp_sched_static;
	} else if (schedule_name == "dynamic") {
		sched = omp_sched_dynamic;
	} else if (schedule_name == "guided") {
		sched = omp_sched_guided;
	} else if (schedule_name == "auto") {
		sched = omp_sched_guided;
	} else {
		throw schedule_name;
	}
	omp_set_schedule(sched, chunk_size);
}

template<typename T>
bool contains(std::vector<T> &vec, T const& value) {
	return std::find(vec.begin(), vec.end(), value) != vec.end();
}

static uint32_t grain_to_color(uint32_t grain) {
	if (grain == 0) return 0x000000;
	if ((grain & 0x3f) == 0) return 0xffffff;
	uint32_t r = ((grain >> 0) << 6 & 0xc0) | 0x3f;
	uint32_t g = ((grain >> 2) << 6 & 0xc0) | 0x3f;
	uint32_t b = ((grain >> 4) << 6 & 0xc0) | 0x3f;
	return r << 16 | g << 8 | b << 0;
}

template<typename F = double, typename M = Material<F, S>>
void compress_scale(ArgumentParser &ap, M &material, F poisson, F T, F dt) {
	F sx = static_cast<F>(1.0)+ap.get<F>("compress-x");
	F sy = static_cast<F>(1.0)+ap.get<F>("compress-y");
	for (auto& c : material.realCells()) {
		c.x = c.x * sx;
		c.y = c.y * sy;
	}
	material.calculateMinMax();
}

template<typename F = double, typename M = Material<F, S>>
void compress_strain(ArgumentParser &ap, M &material, F poisson, F T, F dt) {
	F strain_rate = ap.get<F>("strain-rate");
	F strain = strain_rate * T;
	F i_exp_strain = static_cast<F>(1.0)/exp(strain);

	//F y_size_def = y_size / exp(strain);
	//F x_size_def = (x_size * y_size) / y_size_def;
	for (auto& c : material.realCells()) {
		c.y = c.y * i_exp_strain ;
		c.x = (c.x * c.y) / c.y;
	}
	material.calculateMinMax();
}

template<typename F = double, typename M = Material<F, S>>
void makeImage(ArgumentParser &ap, M &material, int i) {
	auto grain_to_color = [](S grain) -> uint32_t {
		if (grain == 0) return 0x000000;
		if ((grain & 0x3f) == 0) return 0xffffff;
		uint32_t r = ((grain >> 0) << 6 & 0xc0) | 0x3f;
		uint32_t g = ((grain >> 2) << 6 & 0xc0) | 0x3f;
		uint32_t b = ((grain >> 4) << 6 & 0xc0) | 0x3f;
		return r << 16 | g << 8 | b << 0;
	};
	std::stringstream ss;
	ss << "compress/" << std::setw(3) << std::setfill('0') << i << ".bmp";
	int width = 1000;
	std::filesystem::create_directories("compress");
	std::cout << ss.str() << std::endl;
	writeDotsBMP2D<F, rca::Grain, M>(material, ss.str(), width, width, grain_to_color, 0.0, 0.0, width/2, 0, false);
}

template<typename F = double, typename S, typename M>
void correctnessImage(ArgumentParser &ap, M &material) {
	auto path = ap.get<std::string>("correctness-image", "");
	if (path == "") {
		return;
	}
	auto grain_to_color = [](S grain) -> uint32_t {
		if (grain == 0) return 0x000000;
		if ((grain & 0x3f) == 0) return 0xffffff;
		uint32_t r = ((grain >> 0) << 6 & 0xc0) | 0x3f;
		uint32_t g = ((grain >> 2) << 6 & 0xc0) | 0x3f;
		uint32_t b = ((grain >> 4) << 6 & 0xc0) | 0x3f;
		return r << 16 | g << 8 | b << 0;
	};
	int width = 512;
	writeDotsBMP2D<F, S, M>(material, path, width, width, grain_to_color, 0.0, 0.0, width/2, 0, false);
}

std::string papi_results_as_json(unordered_map<std::string, long long> const& results) {
	std::stringstream ss;
	ss << "{";
	size_t i = 0;
	for (auto [event_name, value] : results) {
		ss << "\"" << event_name << "\":\"" << value << "\"";
		if (i+1 < results.size()) {
			ss << ",";
		}
		i++;
	}
	ss << "}";
	return ss.str();
}

void fill_cache(long bytes) {
	if (bytes == 0) {
		return;
	}
	size_t doubles = 1+bytes/sizeof(double);
	vector<double> v(doubles, 0.0);
	for (auto& d : v) {
		d = uniform_random<double>(0.0, 1.0);
	}
	double S = 0.0;
	for (auto& d : v) {
		S += d;
	}
	S /= static_cast<double>(v.size());
	std::cerr << "filled cache (" << S << ")" << flush << std::endl;
}

template<typename F = double>
std::shared_ptr<Neighbourhood<F>> construct_neighbourhood(ArgumentParser &ap, std::string &csv_head, std::string &csv_value) {
	auto nt = ap.get<std::string>("neighbourhood");
	auto dims = ap.get<long>("dimensions", 2);
	std::stringstream ss;
	if (nt == "circle") {
		auto radius = ap.get<F>("radius");
		ss << radius;
		csv_head = "radius";
		csv_value = ss.str();
		if (dims == 2) {
			return std::make_shared<Ellipse<F>>(radius, radius);
		} else if (dims == 3) {
			return std::make_shared<Ellipse<F>>(radius, radius, radius);
		} else {
			throw "bad dimensions";
		}
	} else if (nt == "ellipse") {
		auto ra = ap.get<F>("ra");
		auto rb = ap.get<F>("rb");
		if (dims == 2) {
			csv_head = "ra,rb";
			ss << ra << "," << rb;
			csv_value = ss.str();
			return std::make_shared<Ellipse<F>>(ra, rb);
		} else if (dims == 3) {
			auto rc = ap.get<F>("rc");
			csv_head = "ra,rb,rc";
			ss << ra << "," << rb << "," << rc;
			csv_value = ss.str();
			return std::make_shared<Ellipse<F>>(ra, rb, rc);
		} else {
			throw "bad dimensions";
		}
		} else if (nt == "random-ellipse") {
		auto ra = ap.get<F>("ra");
		auto rb = ap.get<F>("rb");
		if (dims == 2) {
			csv_head = "ra,rb";
			ss << ra << "," << rb;
			csv_value = ss.str();
			return std::make_shared<EllipseRandomRotation<F>>(ra, rb);
		} else {
			throw "bad dimensions";
		}
	} else {
		throw "no such neighbourhood";
	}
}

template<typename F = double>
Extents<F> construct_extents(ArgumentParser &ap) {
	auto dims = ap.get<long>("dimensions", 2);
	auto x_min = ap.get<F>("x-min");
	auto x_max = ap.get<F>("x-max");
	auto y_min = ap.get<F>("y-min");
	auto y_max = ap.get<F>("y-max");
	if (dims == 2) {
		Extents<F> extents(x_min, x_max, y_min, y_max);
		return extents;
	} else if (dims == 3) {
		auto z_min = ap.get<F>("z-min");
		auto z_max = ap.get<F>("z-max");
		Extents<F> extents(x_min, x_max, y_min, y_max, z_min, z_max);
		return extents;
	} else {
		throw "dimensions unrecognized";
	}
}

template<typename F = double>
std::shared_ptr<DistributionGenerator<F>> construct_cell_generator(ArgumentParser &ap, Extents<F> &extents) {
	auto dims = ap.get<long>("dimensions", 2);
	auto cell_generator = ap.get<std::string>("cell-generator");
	if (dims == 2) {
		auto cell_count = ap.get<long>("cell-count");
		if (cell_generator == "grid") {
			auto n = static_cast<long>(sqrt(static_cast<F>(cell_count)));
			return std::make_shared<GridGenerator<F>>(extents, n, n);
		} else if (cell_generator == "random") {
			return std::make_shared<UniformGenerator<F>>(cell_count, extents);
		}else if (cell_generator == "constant") {
			return std::make_shared<ConstantGenerator<F>>(cell_count, extents);
		} else {
			throw "no such cell generator";
		}
	} else if (dims == 3) {
		auto cell_count = ap.get<long>("cell-count");
		if (cell_generator == "grid") {
			auto n = static_cast<long>(pow(static_cast<F>(cell_count), 1.0/3.0));
			return std::make_shared<GridGenerator<F>>(extents, n, n, n);
		} else if (cell_generator == "random") {
			return std::make_shared<UniformGenerator<F>>(cell_count, extents);
		}else if (cell_generator == "constant") {
			return std::make_shared<ConstantGenerator<F>>(cell_count, extents);
		} else {
			throw "no such cell generator";
		}
	} else {
		throw "dimensions unrecognized";
	}
}

template<typename F = double, typename M = Material<F, S>>
int selectGrains(M& material, ArgumentParser &ap, unsigned seed = 0) {
		auto grain_arrangement = ap.get<std::string>("grain-arrangement", "random");
		if (grain_arrangement == "random") {
			auto grain_count = ap.get<int>("grain-count", 100);
			material.changeRandomStates(grain_count, [](unsigned i) -> std::tuple<Group, S> { return std::make_tuple<Group, S>(1, i+1); }, seed);
			return grain_count;
		} else if (grain_arrangement == "grid") {
			auto grain_columns = ap.get<int>("grain-grid-columns", 10);
			auto grain_rows = ap.get<int>("grain-grid-rows", 10);
			auto extents = material.getExtents();

			F dx = extents.width() / static_cast<F>(grain_columns);
			F dy = extents.height() / static_cast<F>(grain_rows);
			F Y = extents.y_min + static_cast<F>(0.5)*dy;
			int i = 1;
			for (int y = 0; y < grain_rows; y++, Y += dy) {
				F X = extents.x_min + static_cast<F>(0.5)*dx;
				for (int x = 0; x < grain_columns; x++, X += dx) {
					selectCellClosestTo<F, S, M>(material, Grain(i), X, Y);
					i++;
				}
			}
			return grain_columns*grain_rows;
		} else {
			throw "unknown grain arrangement";
		}
}

template<typename F, typename M>
int benchmark_mpi(ArgumentParser &ap) {
	Timer<> work, whole;
	MPIProcess<F>::init(&ap.argc, &ap.argv);
	auto omp_num_threads = ap.get<unsigned>("omp-num-threads", 1);
	omp_set_num_threads(omp_num_threads);
	setOmpSchedule(ap);
	auto globalExtents = construct_extents<F>(ap);
	std::tuple<int, int> mpi_grid = {
		ap.get<int>("mpi-grid-width", 0),
		ap.get<int>("mpi-grid-height", 0),
	};
	auto process = MPIProcess<F>::configure(globalExtents, mpi_grid);
	auto dims = ap.get<int>("dimensions", 2);
	std::string a = "a";
	auto neighbourhood = construct_neighbourhood<F>(ap, a, a);
	auto output = ap.get<std::string>("output");
	auto rule = ap.get<string>("rule");
	auto fill_cache_bytes = ap.get<long>("fill-cache", 0);
	auto exchange_edges = ap.get<int>("mpi-exchange-edges", 1);
	size_t changes = 1;
	Timer<> t;
	unsigned image_width, image_height;
	std::string image_directory;
	vector<double> preps;
	vector<double> steps;
	vector<double> ex;
	preps.reserve(1000);
	steps.reserve(1000);
	ex.reserve(1000);
	if (rule != "grain-growth") {
		throw "rule not supported";
	}

	whole.start();
	M material(process.extents);
	auto dg = construct_cell_generator<F>(ap, process.extents);
	material.setNeighbourhoodAndDimensions(neighbourhood, dims);
	material.generate_using(dg.get());
	auto grain_count = selectGrains<F, M>(material, ap, process.id);
	material.ruleFactory = std::make_shared<rca::gg::GrainGrowthRuleFactory<F, S>>();
	std::unique_ptr<Simulation<F, S, M>> simulation = get_simulation<F, M>(ap, material);
	auto cell_count = material.numAllCells();

	fill_cache(fill_cache_bytes);
	process.barrier("setup");

	auto do_next_step = [&]() -> bool {
		bool e = false;
		t.start();
		if (exchange_edges) {
			e = material.exchangeEdges(process, changes);
		}
		t.stop();
		ex.push_back(t.time<double>());
		if (!exchange_edges) {
			e = changes > 0;
		}
		return e;
	};

	auto write_image = [&](unsigned step) {
		auto scale = image_width/std::max(globalExtents.width(), globalExtents.height());
		auto center_x = (globalExtents.x_min + globalExtents.x_max)*0.5;
		auto center_y = (globalExtents.y_min + globalExtents.y_max)*0.5;
		std::stringstream ss;
		ss << image_directory << "/step-" << process.id << "-" << setw(3) << setfill('0') << step << ".bmp";
		std::cout << ss.str() << std::endl;
		writeDotsBMP2D<F, S, M>(material, "final", image_width, image_height,
				grain_to_color, center_x, center_y, scale, 0, false);
	};
	unsigned step = 0;

	if (output == "images") {
		image_width = ap.get<unsigned>("image-width");
		image_height = ap.get<unsigned>("image-height");
		image_directory = ap.get<std::string>("images-directory");
		std::filesystem::create_directories(image_directory);
		//write_image(0);
	}
	work.start();
	while (do_next_step()) {
		t.start();
		simulation->prepare();
		t.stop();
		preps.push_back(t.time<double>());

		t.start();
		changes = simulation->step();
		t.stop();
		steps.push_back(t.time<double>());
		step++;
		// if (output == "images") {
		// 	//write_image(0);
		// }
	}
	if (output == "images") {
		image_width = ap.get<unsigned>("image-width");
		image_height = ap.get<unsigned>("image-height");
		image_directory = ap.get<std::string>("images-directory");
		std::filesystem::create_directories(image_directory);
		write_image(0);
		std::cout << "gsdgdfgsdgdfggsdgsg";
	}
			std::cout << "gsdgdfgsdgdfggsdgsg2";
	work.stop();
	whole.stop();
	if (output == "steps") {
		auto omp_nt = omp_get_num_threads();
		if (process.id == 0) {
			std::cout << "process-id,i,w,h,d,cell-count,grain-count,omp-num-threads,exchange,prep,step" << std::endl << std::flush;
		}
		process.barrier("steps");
		for (int p = 0; p < process.numProcesses; p++) {
			if (p == process.id) {
				for (size_t i = 0; i < preps.size(); i++) {
					std::stringstream ss;
					ss << process.id
						<< "," << (i+1)
						<< "," << process.extents.width()
						<< "," << process.extents.height()
						<< "," << process.extents.depth()
						<< "," << cell_count
						<< "," << grain_count
						<< "," << omp_nt
						<< "," << ex[i]
						<< "," << preps[i]
						<< "," << steps[i]
						<< std::endl;
					std::cout << ss.str() << std::flush;
				}
			}
			process.barrier("step");
		}
	}
	process.barrier("stop");
	MPIProcess<F>::unconfigure();
	if (output == "work-time") {
		if (process.id == 0) {
			std::cout << work.time<double>() << std::endl << std::flush;
		}
	} else if (output == "whole-time") {
		if (process.id == 0) {
			std::cout << whole.time<double>() << std::endl << std::flush;
		}
	}
			std::cout << "gsdgdfgsdgdfggsdgsg3";
	return 0;
}

template<typename F, typename M>
int benchmark(ArgumentParser &ap) {
	Timer<> t;
	Papi papi;
	std::string neigh_csv_head, neigh_csv_value;
	vector<double> preps;
	vector<double> steps;
	preps.reserve(1000);
	steps.reserve(1000);

	auto neighbourhood = construct_neighbourhood<F>(ap, neigh_csv_head, neigh_csv_value);
	auto extents = construct_extents<F>(ap);
	auto dg = construct_cell_generator<F>(ap, extents);

	auto rule = ap.get<string>("rule");
	auto output = ap.get<string>("output");
	auto method = ap.get<std::string>("method");
	auto fill_cache_bytes = ap.get<long>("fill-cache", 0);
	auto dims = ap.get<long>("dimensions", 2);

	auto variant = ap.get<std::string>("variant");

	M material(extents);
	material.setNeighbourhoodAndDimensions(neighbourhood, dims);
	material.generate_using(dg.get());
	auto cell_count = material.numAllCells();

	if (rule == "grain-growth") {
		auto grain_count = selectGrains<F, M>(material, ap);

		if ( variant == "most") {
			material.ruleFactory = std::make_shared<gg::GrainGrowthRuleFactory<F, S>>(gg::Variant::Most);
			fill_cache(fill_cache_bytes);
		}else if(variant == "first"){
			std::cout<<"first         fffffffiiiirsttttt2";
			material.ruleFactory = std::make_shared<gg::GrainGrowthRuleFactory<F, S>>(gg::Variant::First);
		}else {
			material.ruleFactory = std::make_shared<gg::GrainGrowthRuleFactory<F, S>>();
		}

		std::unique_ptr<Simulation<F, S, M>> simulation = get_simulation<F, M>(ap, material);
		size_t changes = 0;

		if (output == "steps") {
			do {
				t.start();
				simulation->prepare();
				t.stop();
				preps.push_back(t.time<double>());

				t.start();
				changes = simulation->step();
				t.stop();
				steps.push_back(t.time<double>());
			} while (changes > 0);
			std::stringstream ss_head;
			ss_head << "i,w,h,d,cell-count,grain-count," << neigh_csv_head << ",method,prep,step" << std::endl;
			std::cout << ss_head.str();
			for (size_t i = 0; i < preps.size(); i++) {
				std::stringstream ss;
				ss << (i+1)
					<< "," << extents.width()
					<< "," << extents.height()
					<< "," << extents.depth()
					<< "," << cell_count
					<< "," << grain_count
					<< "," << neigh_csv_value
					<< "," << method
					<< "," << preps[i]
					<< "," << steps[i] << endl;
				std::cout << ss.str();
			}
		} else if (output == "images") {
			auto width = ap.get<unsigned>("image-width");
			auto height = ap.get<unsigned>("image-height");
			auto directory = ap.get<std::string>("images-directory");
			std::filesystem::create_directories(directory);

			auto write_image = [&](unsigned step) {
				auto scale = width/extents.width();
				auto center_x = (extents.x_min + extents.x_max)*0.5;
				auto center_y = (extents.y_min + extents.y_max)*0.5;
				std::stringstream ss;
				ss << directory << "/final.bmp";
				std::cout << ss.str() << std::endl;
				writeDotsBMP2D<F, S, M>(material, ss.str(), width, height,
						grain_to_color, center_x, center_y, scale, 0, false);
			};
			unsigned step = 0;
			write_image(0);
			do {
				simulation->prepare();
				changes = simulation->step();
				write_image(++step);
			} while (changes > 0);
		} else if (output == "papi-steps") {
			auto papi_events = ap.get<std::string>("papi-events");
			Papi papi;
			papi.init();
			for (auto event : split(papi_events, "+")) {
				papi.addEvent(event.c_str());
			}
			vector<unordered_map<std::string, long long>> prep_papis;
			vector<unordered_map<std::string, long long>> step_papis;
			prep_papis.reserve(1000);
			step_papis.reserve(1000);
			do {
				t.start();
				papi.reset();
				papi.start();
				simulation->prepare();
				papi.stop();
				t.stop();

				prep_papis.push_back(papi.results());
				preps.push_back(t.time<double>());

				t.start();
				papi.reset();
				papi.start();
				changes = simulation->step();
				papi.stop();
				t.stop();

				step_papis.push_back(papi.results());
				steps.push_back(t.time<double>());
			} while (changes > 0);

			std::stringstream ss_head;
			ss_head << "i,w,h,d,cell-count,grain-count," << neigh_csv_head << ",method,prep,step";
			for (auto event : split(papi_events, "+")) {
				ss_head << "," << "prep:" << event;
			}
			for (auto event : split(papi_events, "+")) {
				ss_head << "," << "step:" << event;
			}
			ss_head << std::endl;
			std::cout << ss_head.str();
			for (size_t i = 0; i < preps.size(); i++) {
				std::stringstream ss;
				ss << (i+1)
					<< "," << extents.width()
					<< "," << extents.height()
					<< "," << extents.depth()
					<< "," << cell_count
					<< "," << grain_count
					<< "," << neigh_csv_value
					<< "," << method;
				ss << "," << preps[i] << "," << steps[i];
				for (auto event : split(papi_events, "+")) {
					ss << "," << prep_papis[i][std::string(event)];
				}
				for (auto event : split(papi_events, "+")) {
					ss << "," << step_papis[i][std::string(event)];
				}
				ss << std::endl;
				std::cout << ss.str() << std::flush;
			}
		} else {
			throw "no such output";
		}
	} else if (rule == "gg-compression") {
		F T = 0.0;
		F dt = 0.01;
		F poisson = 0.3;
		auto grain_count = selectGrains<F, M>(material, ap);

		
		if ( variant == "most") {
			material.ruleFactory = std::make_shared<gg::GrainGrowthRuleFactory<F, S>>(gg::Variant::Most);
			fill_cache(fill_cache_bytes);
		}else if(variant == "first"){
			std::cout<<"first         fffffffiiiirsttttt";
			material.ruleFactory = std::make_shared<gg::GrainGrowthRuleFactory<F, S>>(gg::Variant::First);
		}else {
			material.ruleFactory = std::make_shared<gg::GrainGrowthRuleFactory<F, S>>();
		}

		fill_cache(fill_cache_bytes);

		auto method = ap.get<std::string>("method");
		unique_ptr<Simulation<F, S, M>> simulation = get_simulation<F, M>(ap, material);
		size_t changes = 0;
		int i = 0;
		std::stringstream ss;
		if (output == "steps") {
			ss << "i,w,h,d,t,cell-count,grain-count," << neigh_csv_head << ",method,prep,step" << std::endl;
			std::cout << ss.str();
		} else if (output == "images") {
			makeImage<F, M>(ap, material, i);
		}
		auto compress_method = ap.get<std::string>("compress-method");
		do {
			if (compress_method == "scale") {
				compress_scale(ap, material, poisson, T, dt);
			} else if (compress_method == "strain") {
				compress_strain(ap, material, poisson, T, dt);
			} else {
				throw "no such compress method";
			}
			T += dt;

			t.start();
			simulation->prepare();
			t.stop();
			auto prep = t.time<double>();

			t.start();
			changes = simulation->step();
			t.stop();
			auto step = t.time<double>();

			i++;
			if (output == "steps") {
				std::stringstream ss;
				ss << i
					<< "," << extents.width()
					<< "," << extents.height()
					<< "," << extents.depth()
					<< "," << T
					<< "," << cell_count
					<< "," << grain_count
					<< "," << neigh_csv_value
					<< "," << method;
				ss << "," << prep << "," << step << std::endl;
				std::cout << ss.str() << std::flush;
			} else if (output == "images") {
				makeImage<F, M>(ap, material, i);
			}
		} while (changes > 0);
	} // rule
	correctnessImage<F, S, M>(ap, material);
	return 0;
}

void display_usage() {
	std::cout
		<< "usage: benchmark [options]..." << std::endl
		<< "a single option is in form: name=value" << std::endl
		<< "unrecognized options are ignored" << std::endl
		<< "available options:" << std::endl
		<< "    floating-precision=int" << std::endl
		<< "        0 = float" << std::endl
		<< "        1 = use double (default)" << std::endl
		<< "        2 = use long double" << std::endl
		<< "    material-version=string  material version to use" << std::endl
		<< "        mpi                  dedicated for MPI" << std::endl
		<< "        state-changed-flag (latest, default)" << std::endl
		<< "        index-not-id" << std::endl
		<< "        separate-arrays" << std::endl
		<< "        xyz-alone (oldest)" << std::endl
		<< "    rule=string              rule (or simulation type)" << std::endl
		<< "        grain-growth" << std::endl
		<< "        gg-compression" << std::endl
		<< "    output=string            what to output" << std::endl
		<< "        steps                output CSV with stats of each step" << std::endl
		<< "        papi-steps           output CSV with stats of each step" << std::endl
		<< "        images               output BMP file for each step" << std::endl
		<< "    neighbourhood=string     type of neighbourhood" << std::endl
		<< "        circle" << std::endl
		<< "        ellipse" << std::endl
		<< "        random-ellipse" << std::endl		
		<< "    cell-count=int           number of cells to generate" << std::endl
		<< "    cell-generator=string" << std::endl
		<< "        grid" << std::endl
		<< "    x-min=float              extents of simulation space" << std::endl
		<< "    x-max=float" << std::endl
		<< "    y-min=float" << std::endl
		<< "    y-max=float" << std::endl
		<< "    z-min=float              required for dimensions=3" << std::endl
		<< "    z-max=float              required for dimensions=3" << std::endl
		<< "    dimensions=int" << std::endl
		<< "        2 = 2D" << std::endl
		<< "        3 = 3D" << std::endl
		<< "    periodic=int" << std::endl
		<< "        0 = off (default)" << std::endl
		<< "        1 = on" << std::endl
		<< "    fill-cache=int           overwrite cache using specified number of bytes" << std::endl
		<< "    correctness-image=str    output 512x512 image of final microstructure to specified path - useful for checking if result is correct" << std::endl
		<< "    method=string" << std::endl
		<< "        Basic" << std::endl
		<< "        SortDimension" << std::endl
		<< "            sort-by-dimensions=string" << std::endl
		<< "                x = sort by X (default)" << std::endl
		<< "                y = sort by Y" << std::endl
		<< "                z = sort by Z" << std::endl
		<< "        QuadTree" << std::endl
		<< "            leaf-limit=int   limit for cells within subtree" << std::endl
		<< "        QuadTree2" << std::endl
		<< "            leaf-limit=int   limit for cells within subtree" << std::endl
		<< "        FixedGrid" << std::endl
		<< "        FixedGridLL" << std::endl
		<< "        FixedGridDirect" << std::endl
		<< "            space-redundancy=float" << std::endl
		<< "                number of cells to allocate in each buckets" << std::endl
		<< "                0.0 = allocate average" << std::endl
		<< "                0.2 = allocate 120% of average" << std::endl
		<< "                -0.2 = allocate 80% of average" << std::endl
		<< "        FixedSubgrid" << std::endl
		<< "            space-redundancy=float" << std::endl
		<< "            subgrid-radius=int" << std::endl
		<< "                \"radius\" of grid neighbourhood. 1 is default" << std::endl
		<< "        FixedGridMutex" << std::endl
		<< "            space-redundancy=float" << std::endl
		<< "        FixedGridGroup" << std::endl
		<< "            space-redundancy=float" << std::endl
		<< "        FixedGridGroup2" << std::endl
		<< "            space-redundancy=float" << std::endl
		<< "        FixedGridGroup3" << std::endl
		<< "            space-redundancy=float" << std::endl
		<< "        FixedGridGroup2TH" << std::endl
		<< "            subgrid-radius=int" << std::endl
		<< "    mpi=int" << std::endl
		<< "        0 = run without MPI (default)" << std::endl
		<< "        1 = run with MPI" << std::endl
		<< "            please note that for mpi=1:" << std::endl
		<< "            - this executable must be ran using mpirun/mpiexec" << std::endl
		<< "            - MOST options work per-process. Some of the most important options affected by this are:" << std::endl
		<< "                  cell-count" << std::endl
		<< "                  grain-count" << std::endl
		<< "                  omp-num-threads" << std::endl
		<< "              cell-count=N for P processes will result in NxP cells total" << std::endl
		<< "            - material extents (x-min,x-max...) are GLOBAL" << std::endl
		<< "              each process receives a fragment of it" << std::endl
		<< "available for mpi=1:" << std::endl
		<< "    mpi-grid-width=int       width of MPI grid" << std::endl
		<< "    mpi-grid-height=int      height of MPI grid" << std::endl
		<< "    mpi-exchange-edges=int   0 = off, 1 = on (default)" << std::endl
		<< "available for output=papi-steps:" << std::endl
		<< "    papi-events=string       PAPI event names separated with + sign" << std::endl
		<< "        To find available PAPI event names run 'papi_avail' on your system" << std::endl
		<< "        Use \"Name\" column" << std::endl
		<< "available for output=images:" << std::endl
		<< "    image-width=int" << std::endl
		<< "    image-height=int" << std::endl
		<< "    images-directory=string" << std::endl
		<< "available for rule=grain-growth:" << std::endl
		<< "    grain-arrangement=string pattern in which to position grains" << std::endl
		<< "        random" << std::endl
		<< "        grid" << std::endl
		<< "available for grain-arrangement=random:" << std::endl
		<< "    grain-count=int          number of initial grains to randomize" << std::endl
		<< "available for grain-arrangement=grid:" << std::endl
		<< "    grain-grid-columns=int   number of grain grid columns" << std::endl
		<< "    grain-grid-rows=int      number of grain grid rows" << std::endl
		<< "available for rule=gg-compression:" << std::endl
		<< "    compress-method=string   compression method to use" << std::endl
		<< "        scale" << std::endl
		<< "        strain" << std::endl
		<< "available for neighbourhood=circle:" << std::endl
		<< "    radius=float             radius of the circle" << std::endl
		<< "available for neighbourhood=ellipse:" << std::endl
		<< "    ra=float                 radius along x axis" << std::endl
		<< "    rb=float                 radius along y axis" << std::endl
		<< "available for compress-method=scale:" << std::endl
		<< "    compress-x=float" << std::endl
		<< "    compress-y=float" << std::endl
		<< "available for compress-method=strain:" << std::endl
		<< "    strain-rate=float" << std::endl
		;
}

int main(int argc, char *argv[]) {
			std::cout << "gsdgdfgsdgdfggsdgsg MAIN";
	try {
		ArgumentParser ap(argc, argv);
		if (ap.empty()) {
			std::cout << "no arguments specified" << std::endl;
			display_usage();
			return 0;
		}
		auto sp = ap.get<int>("floating-precision", 1);
		auto mv = ap.get<std::string>("material-version");
		auto periodic = ap.get<int>("periodic", 0);
		auto omp_num_threads = ap.get<unsigned>("omp-num-threads", 1);
		omp_set_num_threads(omp_num_threads);
		setOmpSchedule(ap);
		auto mpi = ap.get<unsigned>("mpi", 0);
		if (mv == "mpi" || mpi == 1) {
			if (sp == 0) {
				return benchmark_mpi<float, MaterialMPI<float, S>>(ap);
			} else if (sp == 1) {
				return benchmark_mpi<double, MaterialMPI<double, S>>(ap);
			} else if (sp == 2) {
				return benchmark_mpi<long double, MaterialMPI<long double, S>>(ap);
			}
		} else if (mv == "state-changed-flag") {
			if (periodic == 0) {
				if (sp == 0) {
					return benchmark<float, Material<float, S, PeriodicDisabled>>(ap);
				} else if (sp == 1) {
					return benchmark<double, Material<double, S, PeriodicDisabled>>(ap);
				} else if (sp == 2) {
					return benchmark<long double, Material<long double, S, PeriodicDisabled>>(ap);
				}
			} else {
				if (sp == 0) {
					return benchmark<float, Material<float, S, PeriodicEnabled>>(ap);
				} else if (sp == 1) {
					return benchmark<double, Material<double, S, PeriodicEnabled>>(ap);
				} else if (sp == 2) {
					return benchmark<long double, Material<long double, S, PeriodicEnabled>>(ap);
				}
			}
		} else {
			throw "unknown material version";
		}
	} catch (std::string const& err) {
		std::cerr << "error: " << err << endl;
		return 1;
	} catch (const char* err) {
		std::cerr << "error: " << err << endl;
		return 1;
	} catch (int const& i) {
		std::cerr << "error code " << i << endl;
		return i;
	}
}

// vim: noai:ts=4:sw=4
