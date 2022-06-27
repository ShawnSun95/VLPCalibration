// This is a program that calibrate the Lambertian parameters a, m and M.
// by Xiao Sun

#include<fstream>
#include "ceres/ceres.h"
#include "glog/logging.h"
#include<cmath>

using ceres::AutoDiffCostFunction;
using ceres::NumericDiffCostFunction;
using ceres::Covariance;
using ceres::CostFunction;
using ceres::Problem;
using ceres::CauchyLoss;
using ceres::HuberLoss;
using ceres::Solver;
using ceres::Solve;
using ceres::ResidualBlockId;
using namespace std;
using namespace Eigen;

const int obnum = 71;

struct VLCResidual {
	VLCResidual(double RSS, double d, int index)
		: RSS_(RSS), d_(d), index_(index) {}

	template <typename T> bool operator()(const T* const a,
                                        const T* const b,
                                        T* residual) const {
		residual[0] = a[index_] * pow(RSS_, b[index_]) - d_;
		residual[0] = residual[0] * RSS_;
		return true;
	}

	private:
	const double RSS_;
	const double d_;
	const int index_;
};

int main(int argc, char** argv) {
	google::InitGoogleLogging(argv[0]);

	if (argc < 2) {
		std::cerr << "parameter missing!\n";
		return 1;
	}
	double *x=new double[obnum*3];
	double *RSS=new double[obnum*5];
	//double *variation=new double[obnum*5];
	//double *std=new double[obnum*5];

	fstream infile;
	fstream outfile;
	infile.open(argv[1], ios::in);
	outfile.open("./output/parameter.txt", ios::out);
	for(int i=0;i<obnum;i++)
	{
		for(int j=0;j<3;j++)
			infile >> x[i*3+j];
		for(int j=0;j<5;j++)
			infile >> RSS[i*5+j];
		//for(int j=0;j<5;j++)
		//	infile >> variation[i*5+j];
		//for(int j=0;j<5;j++)
		//	infile >> std[i*5+j];
	}
	infile.close();

	double a[5] = {8.4981,7.2491,7.0180,7.4279,7.4600};
	double b[5] = {-0.3200,-0.2842,-0.2784,-0.2927,-0.3082};

	//double a[5] = {2,2,2,1,1};
	//double b[5] = {-0.3,-0.3,-0.2,-0.3,-0.3};
	double LED[15]={0.26,0.66,2.6,
	    4.61,0.67,2.6,
	    4.57,4.25,2.6,
	    0.18,4.25,2.6,
	    2.45,2.49,2.6};
	double d;
	Problem problem;
	int runtime_number_of_residuals = 3;
	for (int i = 0; i < 5*obnum; ++i) {
		d = sqrt(pow(LED[i%5*3] - x[i/5*3], 2) 
			+ pow(LED[i%5*3+1] - x[i/5*3+1], 2) 
			+ pow(LED[i%5*3+2] - x[i/5*3+2], 2));
		std::cout << d << '\t'<<RSS[i] << "\n";
		problem.AddResidualBlock(
		new NumericDiffCostFunction<VLCResidual, ceres::CENTRAL, ceres::DYNAMIC, 5, 5>(
		new VLCResidual(RSS[i], d, i%5),ceres::TAKE_OWNERSHIP,1),
		NULL,
		a,b);
	}

	Solver::Options options;
	//options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;
	//options.dogleg_type = ceres::SUBSPACE_DOGLEG;
	options.max_num_iterations = 100;
	options.linear_solver_type = ceres::DENSE_QR;
	options.minimizer_progress_to_stdout = true;

	Solver::Summary summary;
	ceres::Solve(options, &problem, &summary);
	
	vector<ResidualBlockId> residual_blocks;
	problem.GetResidualBlocks(&residual_blocks);
	/* Calculate residual */
	ceres::Problem::EvaluateOptions EvalOpts;
	EvalOpts.num_threads = 8;
	EvalOpts.apply_loss_function = false;
	EvalOpts.residual_blocks = residual_blocks;

	std::vector<double> Residuals;
	problem.Evaluate(EvalOpts, NULL, &Residuals, NULL, NULL);
	fstream outfile2;
	outfile2.open("./output/residuals.txt", ios::out);
	for(int i = 0;  i < Residuals.size(); i++)
	{
		d = sqrt(pow(LED[i%5*3] - x[i/5*3], 2) 
			+ pow(LED[i%5*3+1] - x[i/5*3+1], 2));
		outfile2 << d << '\t' << Residuals.at(i)/RSS[i] << endl;
	}

	std::cout << summary.FullReport() << "\n";
	for (int i = 0; i < 5; ++i) {
		outfile << a[i] << '\t' << b[i]  << endl;
	}
	std::cout << "time counts:"<< summary.jacobian_evaluation_time_in_seconds << endl;
	outfile.close();
	return 0;
}
