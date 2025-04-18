#pragma once

#include "Core.h"
#include <random>
#include <algorithm>

class Sampler
{
public:
	virtual float next() = 0;
};

class MTRandom : public Sampler
{
public:
	std::mt19937 generator;
	std::uniform_real_distribution<float> dist;
	MTRandom(unsigned int seed = 1) : dist(0.0f, 1.0f)
	{
		generator.seed(seed);
	}
	float next()
	{
		return dist(generator);
	}
};

// Note all of these distributions assume z-up coordinate system
class SamplingDistributions
{
public:
	static Vec3 uniformSampleHemisphere(float r1, float r2) {
		float phi = 2.0f * M_PI * r1;         
		float cosTheta = r2;                  
		float sinTheta = sqrt(1.0f - cosTheta * cosTheta);  
		return Vec3(sinTheta * cos(phi), sinTheta * sin(phi), cosTheta);
	}
	static float uniformHemispherePDF(const Vec3 wi) {
		if (wi.z > 0.0f)
			return 1.0f / (2.0f * M_PI);
		else
			return 0.0f;
	}
	static Vec3 cosineSampleHemisphere(float r1, float r2)
	{
		float phi = 2.0f * M_PI * r1;
		float x = cos(phi) * sqrt(r2);
		float y = sin(phi) * sqrt(r2);
		float z = sqrt(1.0f - r2);
		return Vec3(x, y, z);
	}
	static float cosineHemispherePDF(const Vec3 wi)
	{
		if (wi.z > 0.0f)
			return wi.z / M_PI;
		else
			return 0.0f;
	}
	static Vec3 uniformSampleSphere(float r1, float r2)
	{
		float phi = 2.0f * M_PI * r1;      
		float cosTheta = sqrt(r2);
		float sinTheta = sqrt(1.0f - cosTheta * cosTheta);  
		return Vec3(sinTheta * cos(phi), sinTheta * sin(phi), cosTheta);
	}
	static float uniformSpherePDF(const Vec3& wi)
	{
		return 1.0f / (4.0f * M_PI);
	}
	static Vec3 sampleGGX(const Vec3& wo, float alpha, float u1, float u2)
	{
		Vec3 V = Vec3(alpha * wo.x, alpha * wo.y, wo.z).normalize();
		float lensq = V.x * V.x + V.y * V.y;
		Vec3 T1 = lensq > 0.0f ? Vec3(-V.y, V.x, 0).normalize() : Vec3(1, 0, 0);
		Vec3 T2 = V.cross(T1);
		float r = sqrtf(u1);
		float phi = 2.0f * M_PI * u2;
		float t1 = r * cosf(phi);
		float t2 = r * sinf(phi);
		float s = 0.5f * (1.0f + V.z);
		t2 = (1.0f - s) * sqrtf(1.0f - t1 * t1) + s * t2;
		Vec3 Nh = T1 * t1 + T2 * t2 + V * sqrtf(std::max(0.0f, 1.0f - t1 * t1 - t2 * t2));
		Vec3 h = Vec3(alpha * Nh.x, alpha * Nh.y, std::max(0.0f, Nh.z)).normalize();
		return h;
	}
};