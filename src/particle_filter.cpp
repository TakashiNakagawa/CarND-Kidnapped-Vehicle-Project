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
#include <random>
#include <array>
#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[])
{
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1.
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles = 100;
	default_random_engine gen;
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);
	for (auto i = 0; i < num_particles; ++i)
	{
		Particle p;
		p.id = i;
		p.x = dist_x(gen);
		p.y = dist_y(gen);
		p.theta = dist_theta(gen);
		p.weight = 1;
		weights.push_back(1);
		particles.push_back(p);
	}
	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate)
{
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	default_random_engine gen;
	for (auto &p : particles)
	{
		const auto delta_theta = yaw_rate * delta_t;
		if (std::fabs(yaw_rate) < 0.001)
		{
			p.x += velocity * std::cos(p.theta) * delta_t;
			p.y += velocity * std::sin(p.theta) * delta_t;
		}
		else
		{
			p.x += velocity * (std::sin(p.theta + delta_theta) - std::sin(p.theta)) / yaw_rate;
			p.y += velocity * (std::cos(p.theta) - std::cos(p.theta + delta_theta)) / yaw_rate;
			p.theta += delta_theta;
		}
		normal_distribution<double> dist_x(p.x, std_pos[0]);
		normal_distribution<double> dist_y(p.y, std_pos[1]);
		normal_distribution<double> dist_theta(p.theta, std_pos[2]);
		p.x = dist_x(gen);
		p.y = dist_y(gen);
		p.theta = fmod(dist_theta(gen), 2 * M_PI);
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs> &observations)
{
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
	//   implement this method and use it as a helper during the updateWeights phase.
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
								   const std::vector<LandmarkObs> &observations, const Map &map_landmarks)
{
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
	double sum = 0;
	weights.clear();
	auto gauss_norm = (1 / (2 * M_PI * std_landmark[0] * std_landmark[1]));
	std::array<double, 2> denominator =
		{{std::pow(2 * std_landmark[0], 2),
		  std::pow(2 * std_landmark[1], 2)}};
	for (auto &p : particles)
	{
		double weight = 1;
		const auto p_theta_cos = std::cos(p.theta);
		const auto p_theta_sin = std::sin(p.theta);
		std::vector<int> associations;
		std::vector<double> sense_x, sense_y;
		for (auto &ob : observations)
		{
			auto obx = p.x + p_theta_cos * ob.x - p_theta_sin * ob.y;
			auto oby = p.y + p_theta_sin * ob.x + p_theta_cos * ob.y;
			double dist_min = std::numeric_limits<double>::max();
			std::pair<float, float> position;
			const Map::single_landmark_s *nearest_map = nullptr;
			for (const auto &m : map_landmarks.landmark_list)
			{
				if (m.x_f - obx > sensor_range*2 || m.y_f - oby - oby > sensor_range*2)
				{
					continue;
				}
				auto dist = std::pow(m.x_f - obx, 2) + std::pow(m.y_f - oby, 2);
				if (!nearest_map || dist_min > dist)
				{
					dist_min = dist;
					position = std::make_pair(m.x_f, m.y_f);
					nearest_map = &m;
				}
			}
      if(!nearest_map)
      {
        std::cout << "nearest_map null" << std::endl;
        continue;
      }
			auto x_map = position.first;
			auto y_map = position.second;
			associations.push_back(nearest_map->id_i);
			sense_x.push_back(obx);
			sense_y.push_back(oby);

			auto exponent = std::pow(obx - x_map, 2) / denominator[0] + std::pow(oby - y_map, 2) / denominator[1];
			auto ret = gauss_norm * std::exp(-exponent);
			weight *= ret;
		}
    if(associations.empty())
    {
      std::cout << "associations empty" << std::endl;
      continue;
    }
      
		SetAssociations(p, associations, sense_x, sense_y);
		p.weight = weight;
		sum += weight;
		weights.push_back(p.weight);
	}
	if (sum == 0)
	{
		std::cout << "sum = 0" << std::endl;
		std::cout << "particles size = " << particles.size() << std::endl;
	}
}

void ParticleFilter::resample()
{
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	std::random_device rd;
	std::mt19937 gen(rd());
	std::discrete_distribution<> d(std::begin(weights), std::end(weights));
	std::map<int, int> m;
	std::vector<Particle> particles_new;
	for (auto i = 0; i < weights.size(); ++i)
	{
		particles_new.push_back(particles[d(gen)]);
	}
	particles_new.swap(particles);
}

Particle ParticleFilter::SetAssociations(Particle &particle, const std::vector<int> &associations,
										 const std::vector<double> &sense_x, const std::vector<double> &sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	particle.associations = associations;
	particle.sense_x = sense_x;
	particle.sense_y = sense_y;
	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
	copy(v.begin(), v.end(), ostream_iterator<int>(ss, " "));
	string s = ss.str();
	s = s.substr(0, s.length() - 1); // get rid of the trailing space
	return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
	copy(v.begin(), v.end(), ostream_iterator<float>(ss, " "));
	string s = ss.str();
	s = s.substr(0, s.length() - 1); // get rid of the trailing space
	return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
	copy(v.begin(), v.end(), ostream_iterator<float>(ss, " "));
	string s = ss.str();
	s = s.substr(0, s.length() - 1); // get rid of the trailing space
	return s;
}
