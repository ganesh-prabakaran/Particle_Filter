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

#define EPS 0.00001

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
    
	num_particles = 50;
	
	double std_x = std[0];
	double std_y = std[1];
	double std_theta = std[2];
	weights.resize(num_particles);
	// Creates a normal (Gaussian) distribution for gps noise - x, y and theta
	normal_distribution<double> dist_x(x, std_x);
	normal_distribution<double> dist_y(y, std_y);
	normal_distribution<double> dist_theta(theta, std_theta);
	
	particles.resize(num_particles);
	// Add noise to the initial values from GPS.
	for (unsigned int i = 0; i < num_particles; i++){
		particles[i].id = i;
		particles[i].x = dist_x(gen);
	    particles[i].y = dist_y(gen);
	    particles[i].theta = dist_theta(gen);	 
		particles[i].weight = 1.0;
	}

    is_initialized = true;
	
	

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	
	// Predict next location of the particle
	for (unsigned int i = 0; i < num_particles; i++){
		
		if (fabs(yaw_rate) < 0.0001){
			particles[i].x +=  velocity * delta_t * cos(particles[i].theta);
			particles[i].y +=  velocity * delta_t * sin(particles[i].theta);
			
		}
		else{
			particles[i].x += (velocity / yaw_rate) * (sin(particles[i].theta + yaw_rate * delta_t) - sin(particles[i].theta));
			particles[i].y += (velocity / yaw_rate) * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate * delta_t));
			particles[i].theta += yaw_rate * delta_t;
		}
		// Creates a normal (Gaussian) distribution for sensor noise - x, y and theta
		normal_distribution<double> noise_x(0.0, std_pos[0]);
		normal_distribution<double> noise_y(0.0, std_pos[1]);
		normal_distribution<double> noise_theta(0.0, std_pos[2]);
			
		// Add noise to the predicted values.
		 particles[i].x += noise_x(gen);
		 particles[i].y += noise_y(gen);
		 particles[i].theta += noise_theta(gen);	 
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	
	//For each observation 
	for (unsigned int i = 0; i < observations.size(); i++){
		
		double min_dist = numeric_limits<double>::max();
		int closest_pair = -1;
		double o_x = observations[i].x;
		double o_y = observations[i].y;
		double apr_x, apr_y;
		
		//For each landmark in range
		for (unsigned int j = 0; j < predicted.size(); j++){

			double pr_x = predicted[j].x;
			double pr_y = predicted[j].y;
			
			double cur_dist = dist(o_x, o_y,pr_x, pr_y );
			// Find the closest landmark
			if (cur_dist < min_dist){
				min_dist = cur_dist;
				closest_pair = j;
				apr_x = pr_x;
				apr_y = pr_y;
			}
		}
		// Associate the nearest landmark to the observation.
		observations[i].id = closest_pair;
		
	}

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
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
	
	
	// For each particle
	for (unsigned int n = 0; n < num_particles; n++){
		
		vector<LandmarkObs> inrange_map_landmarks;
		
		double p_x = particles[n].x;
		double p_y = particles[n].y;
		double p_theta = particles[n].theta;
		
		// Filter landmarks in range.
		for (unsigned int i = 0; i < map_landmarks.landmark_list.size(); i++){
			double lm_x = map_landmarks.landmark_list[i].x_f;
			double lm_y = map_landmarks.landmark_list[i].y_f;
			int lm_id = map_landmarks.landmark_list[i].id_i;
			double curr_dist = dist(lm_x, lm_y, p_x, p_y );
			if (curr_dist <= sensor_range){
			    inrange_map_landmarks.push_back(LandmarkObs{lm_id, lm_x, lm_y});
			}
		}
		
		// Transform observations in particle coordinates into map coordinates. 
		vector<LandmarkObs> mapped_observations;
		LandmarkObs t_obs, assoc_lm_obs;
		for(unsigned int i =0; i < observations.size(); i++){
		   double o_x = observations[i].x;
		   double o_y = observations[i].y;
		   t_obs.x = p_x + (cos(p_theta) * o_x) - (sin(p_theta) * o_y);
		   t_obs.y = p_y + (sin(p_theta) * o_x) + (cos(p_theta) * o_y);
		   t_obs.id = observations[i].id;
           mapped_observations.push_back(t_obs);
		}
	    
		// Associate each observation with a nearest landmark
        dataAssociation(inrange_map_landmarks, mapped_observations);
        
		vector<int> associations;
		vector<double> sense_x;
		vector<double> sense_y;
		
		// Assign the associations to the particle
		for (unsigned int i = 0; i < mapped_observations.size(); i++){
			associations.push_back(inrange_map_landmarks[mapped_observations[i].id].id);
			sense_x.push_back(mapped_observations[i].x);
			sense_y.push_back(mapped_observations[i].y);
		}
		
		SetAssociations(particles[n], associations, sense_x, sense_y);
		
		
		// Weight calculation
        double std_x = std_landmark[0];
        double std_y = std_landmark[1];
		   
		double gauss_norm = 1 / (2 * M_PI * std_x * std_y);
		particles[n].weight = 1.0;		
		double weight;
		
		for(unsigned int i = 0; i < mapped_observations.size(); i++){
           
		   
		   int associated_id = mapped_observations[i].id;
		   double mo_x = mapped_observations[i].x;
		   double mo_y = mapped_observations[i].y;
		   double alm_x, alm_y;
		   
		   assoc_lm_obs = inrange_map_landmarks[associated_id];
		   
		   // Calculate weight using multi variate gaussain distribution
		   double x_diff = mo_x - assoc_lm_obs.x;
		   double y_diff = mo_y - assoc_lm_obs.y;
		   double exponent = ((x_diff * x_diff) / (2 * std_x * std_x)) + ((y_diff * y_diff) / (2 * std_y * std_y));
		   
		   weight = gauss_norm * exp(-exponent); 
		   if (weight == 0){
			   particles[n].weight *= EPS;
		   }
		   else{
			   particles[n].weight *= weight; 
		   }
	    }	   
		//Assign the weight to the particle.   
		particles[n].weight = weight; 
		weights[n]=weight;
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
  
  vector<Particle> resampled; 
  
  for (unsigned int i = 0; i < num_particles; ++i){
	  std::discrete_distribution<int> d(weights.begin(), weights.end());
	  resampled.push_back(particles[d(gen)]);
  }
  particles = resampled;
  
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

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
