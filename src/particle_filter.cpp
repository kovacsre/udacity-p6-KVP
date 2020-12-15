/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>
#include <cmath>

#include "helper_functions.h"

using namespace std;
using std::string;
using std::vector;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 10;  // TODO: Set the number of particles


  //create normal distributions for sensor noise
  default_random_engine gen;
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_th(theta, std[2]);

  //initialize particles
  for (int i=0; i<num_particles; i++){
	 Particle p;
	 p.id = i;
	 p.x = dist_x(gen);
	 p.y = dist_y(gen);
	 p.theta = dist_th(gen);
	 p.weight = 1.0;

	 particles.push_back(p);

 }
   is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
    // particles structot kell update-elni, x,y,theta-t a marha nagy kepletbol kell kiszamolni a leckeben
	//kell zaj is a measurementekre (mar van rajtuk, nem?)

	//create normal distributions
	  default_random_engine gen;
	  normal_distribution<double> dist_x(0, std_pos[0]);
	  normal_distribution<double> dist_y(0, std_pos[1]);
	  normal_distribution<double> dist_th(0, std_pos[2]);

	  for (int i =0; i<num_particles; i++){
		  //get the actual values
		  double x0 = particles[i].x;
		  double y0 = particles[i].y;
		  double th0 = particles[i].theta;
		  //update particle - avoid dividing by zero
		  if (fabs(yaw_rate) > 0.0001){
			  particles[i].x = x0 + velocity/yaw_rate * (sin(th0 + yaw_rate * delta_t) - sin(th0));
			  particles[i].y = y0 + velocity/yaw_rate * (cos(th0) - cos(th0 + yaw_rate * delta_t));
		  } else{
			  particles[i].x = x0 + velocity * delta_t * cos(th0);
			  particles[i].y = y0 + velocity * delta_t * sin(th0);
		  }
		  particles[i].theta = th0 + yaw_rate * delta_t;
		  //add random noise
		  particles[i].x += dist_x(gen);
		  particles[i].y += dist_y(gen);
		  particles[i].theta += dist_th(gen);


	  }

	 // cout<< "prediction done" << endl;
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
	//vegig kell menni az observations array-en es mindegyik tavolsagat megnezni a predictedtol
	//amelyik a legkisebb, azt az idt kell eltarolni

	for(unsigned int i=0; i<observations.size(); i++){
		double min_distance = 999;
		//current observation
		LandmarkObs obs = observations[i];

		int lm_id = -1;

		for(unsigned int j=0; j<predicted.size(); j++){
			LandmarkObs pre = predicted[j];
			double distance = dist(pre.x, pre.y, obs.x, obs.y);

			if(fabs(distance) < min_distance){
				min_distance = distance;
				lm_id = pre.id;
			}
		}

		//the observation's id shall be the same as the id of the nearest landmark
		observations[i].id = lm_id;
	}
	//cout<< "data association done" << endl;
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */

	vector<double> w; //vector to store all weights - for debugging purpose

	//to be done for each particle
	for(unsigned int i=0; i<particles.size(); i++){
		//particle coordinates
		double px = particles[i].x;
		double py = particles[i].y;
		double ptheta = particles[i].theta;

		//landmarkobs vector to store landmarks in range
		vector<LandmarkObs> lm_in_range;
		LandmarkObs landm;

		//for each landmark
		for(unsigned int j=0; j<map_landmarks.landmark_list.size(); j++){
			//landmark coordinates
			float lx = map_landmarks.landmark_list[j].x_f;
			float ly = map_landmarks.landmark_list[j].y_f;
			int l_id = map_landmarks.landmark_list[j].id_i;

			//distance of landmarks from particle
			double lm_dist = dist(px, py, lx, ly);

			//if landmark is in range
			if(fabs(lm_dist) <= sensor_range){
				landm.x = lx;
				landm.y = ly;
				landm.id = l_id;
				//store values in vector
				lm_in_range.push_back(landm);
			}
		}

		//transform observations from vehicle coordinates to map coordinates
		vector<LandmarkObs> observations_mapcoord;
		LandmarkObs observation;

		for(unsigned int k=0; k<observations.size(); k++){
			observation.x = px + ((cos(ptheta) * observations[k].x) - (sin(ptheta) * observations[k].y));
			observation.y = py + ((sin(ptheta) * observations[k].x) + (cos(ptheta) * observations[k].y));
			observation.id = observations[k].id;
			//store values in vector
			observations_mapcoord.push_back(observation);
		}

		//data association on landmarks in range and transformed observations - to get nearest landmark to particle
		dataAssociation(lm_in_range, observations_mapcoord);

		double obs_x, obs_y, mu_x, mu_y;

		//get the coordinates of the associated landmark
		for (unsigned int m=0; m<observations_mapcoord.size(); m++){
			obs_x = observations_mapcoord[m].x;
			obs_y = observations_mapcoord[m].y;
			int associated_pred = observations_mapcoord[m].id;

			for(unsigned int n=0; n<lm_in_range.size(); n++){
				if(lm_in_range[n].id == associated_pred){
					mu_x = lm_in_range[n].x;
					mu_y = lm_in_range[n].y;
				}
			}
		}
		//use the multivariate Gaussian distribution to update the weight of each particle
		double weight = multiv_prob(std_landmark[0], std_landmark[1], obs_x, obs_y, mu_x, mu_y);

		particles[i].weight = weight;
		w.push_back(weight);
	}
	//cout<< "weight update done" << endl;
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */


	//create the new particles
	vector<Particle> new_particles;

	//get the weights
	vector<double> weights;
	for(int i=0; i<num_particles; i++){
		weights.push_back(particles[i].weight);
	}


	//resampling wheel rewritten from the python code
	//get the starting index
	default_random_engine gen;
	//uniform distribution
	uniform_int_distribution<int> d(0, num_particles-1);
	int index = d(gen);
	double max_weight = *max_element(weights.begin(), weights.end());
	double beta = 0.0;

	//create a uniform distribution for the random value to be added to beta later
	if(max_weight != 0.0){
		uniform_real_distribution<double> u(0.0, 2*max_weight);

		for(unsigned int j=0; j<particles.size(); j++){
			double offset = u(gen);
			beta += offset;
			while(beta > weights[index]){
				beta -= weights[index];
				index = (index + 1) % num_particles;
			}
			new_particles.push_back(particles[index]);
		}
		particles = new_particles;
	}

	//cout<< "resampling done" <<endl;
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}
