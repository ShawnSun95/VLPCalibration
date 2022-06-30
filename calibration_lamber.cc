// This is a program that calibrate the Lambertian parameters a, m and M.
// by Xiao Sun

#include<fstream>
#include "ceres/ceres.h"
#include "glog/logging.h"
#include<cmath>

using ceres::AutoDiffCostFunction;
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
using namespace ceres;

class VLPCalibration : public SizedCostFunction<1,5,1,5> {
	public:
		VLPCalibration(double RSS, Vector3d p_pd, Vector3d p_led, double alpha, int index)
			: RSS_(RSS), p_pd_(std::move(p_pd)), p_led_(std::move(p_led)), alpha_(alpha), index_(index) {}

		virtual ~VLPCalibration() {}

		virtual bool Evaluate(double const* const* parameters,
                           double* residuals,
                           double** jacobians) const {
			const double *m = parameters[0];
			const double *M = parameters[1];
			const double *a = parameters[2];

			Vector3d e_PD(0, sin(alpha_), cos(alpha_));
			Vector3d e_LED(0, 0, 1);
			Vector3d LOS = p_led_ - p_pd_;
			double cos_theta = LOS.dot(e_LED) / LOS.norm();
			double cos_phi = LOS.dot(e_PD) / LOS.norm();
			double d2 = LOS.norm()*LOS.norm();
			residuals[0] = a[index_]*(m[index_]+1)*pow(cos_theta, m[index_])*pow(cos_phi, M[0])/d2 - RSS_;
			
			if (!jacobians) return true;
			double* jacobian0 = jacobians[0];
			if (!jacobian0) return true;
			double* jacobian1 = jacobians[1];
			if (!jacobian1) return true;
			double* jacobian2 = jacobians[2];
			if (!jacobian2) return true;

			for (int i=0;i<5;i++){
				jacobian0[i] = 0;
				jacobian2[i] = 0;
			}
			
			jacobian0[index_] = a[index_]*(m[index_]+1)*log(cos_theta)*pow(cos_theta, m[index_])*pow(cos_phi, M[0])/d2
				+a[index_]*pow(cos_theta, m[index_])*pow(cos_phi, M[0])/d2;
			jacobian1[0] = a[index_]*log(cos_phi)*pow(cos_theta, m[index_])*pow(cos_phi, M[0])/d2;
			jacobian2[index_] = pow(cos_theta, m[index_])*pow(cos_phi, M[0])/d2;

			return true;
		}
	private:
		const double RSS_;
		const Vector3d p_pd_;
		const Vector3d p_led_;
		const double alpha_;
		const int index_;
};


// struct VLPCalibration {
// 	VLPCalibration(double RSS, Vector3d p_pd, Vector3d p_led, double alpha, int index)
// 		: RSS_(RSS), p_pd_(std::move(p_pd)), p_led_(std::move(p_led)), alpha_(alpha), index_(index) {}

// 	template <typename T> bool operator()(const T* const m,
//                                         const T* const M,
// 										const T* const a,
//                                         T* residual) const {
// 		Vector3d e_PD(0, sin(alpha_), cos(alpha_));
// 		Vector3d e_LED(0, 0, 1);
// 		Vector3d LOS = p_led_ - p_pd_;
// 		T cos_theta = LOS.dot(e_LED) / LOS.norm();
// 		T cos_phi = LOS.dot(e_PD) / LOS.norm();
// 		residual[0] = a[index_]*pow(cos_theta, m[index_])*pow(cos_phi, M[0])/LOS.norm()/LOS.norm() - RSS_;
// 		return true;
// 	}

// 	private:
// 	const double RSS_;
// 	const Vector3d p_pd_;
// 	const Vector3d p_led_;
// 	const double alpha_;
// 	const int index_;
// };

int main(int argc, char** argv) {
	google::InitGoogleLogging(argv[0]);

	if (argc < 5) {
		std::cerr << "parameter missing!\n";
		return 1;
	}

	int n_angle, n_pos, n_led;
	n_angle = atoi(argv[2]);
	n_pos   = atoi(argv[3]);
	n_led   = atoi(argv[4]);
	double *x     = new double[n_angle*n_pos*3];
	double *RSS   = new double[n_angle*n_pos*n_led];
	double *LED   = new double[n_led*3];
	double *alpha = new double[n_angle];
	double *m     = new double[n_led];
	double *a     = new double[n_led];
	double M, PD_Height;

	fstream infile;
	fstream outfile;
	infile.open(argv[1], ios::in);
	outfile.open("./output/parameter_lamber.txt", ios::out);
	for(int i=0;i<n_angle;i++){
		for(int j=0;j<n_pos;j++){
			for(int k=0;k<3;k++)
				infile >> x[i*n_pos*3+j*3+k];
			for(int k=0;k<n_led;k++)
				infile >> RSS[i*n_pos*n_led+j*n_led+k];
		}
	}
	for(int k=0;k<n_led;k++){
		infile >> LED[k*3] >> LED[k*3+1] >> LED[k*3+2];
	}
	for(int i=0;i<n_angle;i++)
		infile >> alpha[i];
	infile >> PD_Height >> m[0] >> M;
	for(int k=0;k<n_led;k++){
		infile >> a[k];
		m[k] = m[0];
	}

	infile.close();
	vector<Vector3i> inds;

	Problem problem;
	for (int i = 0; i < n_angle; ++i) {
		for (int j = 0; j < n_pos; j++){
			for (int k = 0; k < n_led; k++){
				Vector3d p_LED(LED[k*3], LED[k*3+1], LED[k*3+2]);
				Vector3d p_PD(x[i*n_pos*3+j*3], x[i*n_pos*3+j*3+1], x[i*n_pos*3+j*3+2]);
				Vector3d e_PD(0, sin(alpha[i]), cos(alpha[i]));
				Vector3d e_LED(0, 0, 1);
				Vector3d LOS = p_LED - p_PD;
				double cos_theta = LOS.dot(e_LED) / LOS.norm();
				double cos_phi = LOS.dot(e_PD) / LOS.norm();
				if (cos_theta < 0.5 || cos_phi < 0.5) continue;
				
				CostFunction* cost_function =new VLPCalibration(RSS[i*n_pos*n_led+j*n_led+k], p_PD, p_LED, alpha[i], k);
				problem.AddResidualBlock(cost_function,	nullptr, m, &M, a);

				inds.push_back(Vector3i(i,j,k));
			}
		}
	}

	Solver::Options options;
	options.max_num_iterations = 100;
	// options.linear_solver_type = ceres::DENSE_QR;
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
	outfile2.open("./output/residuals_lamber.txt", ios::out);
	for(int i = 0;  i < Residuals.size(); i++)
	{
		Vector3i ind=inds[i];
		outfile2 << ind(0)+1 << '\t' << ind(1)+1 << '\t' << ind(2)+1 << '\t' <<
			RSS[ind(0)*n_pos*n_led+ind(1)*n_led+ind(2)] << '\t' << Residuals.at(i) << 
			'\t' << Residuals.at(i)/RSS[ind(0)*n_pos*n_led+ind(1)*n_led+ind(2)] << endl;
	}

	std::cout << summary.FullReport() << "\n";
	for (int i = 0; i < 5; ++i) {
		outfile << m[i] << '\t' << M << '\t' << a[i]  << endl;
	}
	std::cout << "time counts:"<< summary.jacobian_evaluation_time_in_seconds << endl;
	outfile.close();
	return 0;
}
