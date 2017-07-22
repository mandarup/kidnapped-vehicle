/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h>
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std_pos[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1.
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles = 100;
	std::default_random_engine gen;
	std::normal_distribution<double> N_x(x,std_pos[0]);
	std::normal_distribution<double> N_y(y,std_pos[1]);
	std::normal_distribution<double> N_theta(theta,std_pos[2]);

	for(int i = 0; i < num_particles; i++){
		Particle particle;
		particle.id = i;
		particle.x = N_x(gen);
		particle.y = N_y(gen);
		particle.theta = N_theta(gen);
		particle.weight = 1;

		particles.push_back(particle);
		particles.push_back(1);
	}
	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	std::default_random_engine gen;

	for(int i = 0; i < num_particles; i++){
		double new_x;
		double new_y;
		double new_theta;

		double x0 = particles[i].x;
		double y0 = particles[i].y;
		double theta0 = particles[i].theta;


		if(yaw_rate == 0){
			new_x = x0 + velocity * (sin(theta0) - sin(theta0));
			new_y = y0 + velocity *(-cos(theta0)+ cos(theta0));
			new_theta = theta0 + yaw_rate * delta_t;
		}
		else{
			new_x = x0 + velocity/yaw_rate * (sin(theta0 + yaw_rate * delta_t)
											  - sin(theta0));
			new_y = y0 + velocity/yaw_rate *(-cos(theta0 + yaw_rate * delta_t)
											  + cos(theta0));
			new_theta = theta0 + yaw_rate * delta_t;
		}

		std::normal_distribution<double> N_x(x,std_pos[0]);
		std::normal_distribution<double> N_y(y,std_pos[1]);
		std::normal_distribution<double> N_theta(theta,std_pos[2]);

		particles[i].x = N_x(gen);
		particles[i].y = N_y(gen);
		particles[i].theta = N_theta(gen);

	}



}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
	//   implement this method and use it as a helper during the updateWeights phase.

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html


	for (int j = 0; j < particles.size(); j++){
		p = particles[j];

		std::vector<int> associations;
		std::vector<double> sense_x;
		std::vector<double> sense_y;

		std::vector<LandmarkObs> trans_obs;
		LandmarkObs obs;
		for (int i = 0; i < observations.size(); i++){
			LandmarkObs trans_ob;
			obs = observations[i];

			// transform from map coordinates to vehicle coordinates
			trans_ob.x = p.x + (obs.x * cos(p.theta) - obs.y * sin(p.theta));
			trans_ob.y = p.y + (obs.x * sin(p.theta) + obs.y * cos(p.theta));
			trans_obs.push_back(trans_ob)
		}
		p.weight = 1.0;

		for(int i; i < trans_obs.size(); i++){
			double closest_dist = sensor_range;
			int association = 0;

			for(int k; k < map_landmarks.size(); k++){
				double landmark_x = map_landmarks.landmark_list[k].x_f;
				double landmark_y = map_landmarks.landmark_list[k].y_f;
				double calc_dist = math::sqrt(math::pow(trans_obs[i].x - landmark_x, 2.0)
											 + math::pow(trans_obs[i].y - landmark_y, 2.0));
				if(calc_dist < closest_dist){
					closest_dist = calc_dist;
					association = k;
				}
			}

			if(association != 0){
				double meas_x = trans_obs[i].x;
				double meas_y = trans_obs[i].y;
				double mu_x = map_landmarks.landmark_list[association].x_f;
				double mu_y = map_landmarks.landmark_list[association].y_f;

				long double multiplier = 1/(2 * math::M_PI * std_landmark[0] * std_landmark[1])
										* math::exp(-( gaussianKernel(meas_x, mu_x, std_landmark[0])
													  + gaussianKernel(meas_y, mu_y, std_landmark[1])));
				if (multiplier > 0){
					p.weight *= multiplier;
				}
			}
			associations.push_back(association + 1);
			sense_x.push_back(trans_obs[i].x);
			sense_y.push_back(trans_obs[i].y);
		}


	}


}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	default_random_engine gen;
	discrete_distribution<int> distribution(weights.begin(), weights.end());
	std::vector<Particle> resample_particles;

	for(int i=0; i < num_particles; i++){
		resample_particles.push_back(particles[distribution(gen)]);
	}
	particles = resample_particles;

}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}


double gaussianKernel(double x, double mu, double sigma){
	return pow((x - mu)/sigma, 2.0)/2.0
}
