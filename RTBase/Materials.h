#pragma once

#include "Core.h"
#include "Imaging.h"
#include "Sampling.h"

#pragma warning( disable : 4244)

class BSDF;
class ShadingData
{
public:
	Vec3 x;
	Vec3 wo;
	Vec3 sNormal;
	Vec3 gNormal;
	float tu;
	float tv;
	Frame frame;
	BSDF* bsdf;
	float t;
	ShadingData() {}
	ShadingData(Vec3 _x, Vec3 n)
	{
		x = _x;
		gNormal = n;
		sNormal = n;
		bsdf = NULL;
	}
};

template<typename T>
inline T clamp(T val, T minVal, T maxVal)
{
	return std::max(minVal, std::min(maxVal, val));
}

class ShadingHelper
{
public:
	static Vec3 reflect(const Vec3& input, const Vec3& normal)
	{
		return -input + normal * 2.0f * Dot(input, normal);
	}
	//GGX
	static float Dggx(Vec3 h, float alpha)
	{
		if (h.z <= 0.0f) return 0.0f;

		float alpha2 = alpha * alpha;
		float cosTheta = std::max(h.z, 1e-6f);
		float cosTheta2 = cosTheta * cosTheta;
		float tan2Theta = 1.0f / cosTheta2 - 1.0f;
		tan2Theta = std::max(tan2Theta, 0.0f); 
		float tmp = alpha2 + tan2Theta;
		float denom = (float)M_PI * cosTheta2 * cosTheta2 * (tmp * tmp);
		if (denom < 1e-12f || !std::isfinite(denom)) {
			return 0.0f;
		}

		return alpha2 / denom;
	}
	//lambda
	static float lambdaGGX(Vec3 wi, float alpha)
	{
		float absTanTheta = fabsf(sqrtf(1.0f - wi.z * wi.z) / wi.z);
		if (std::isinf(absTanTheta)) return 0.0f;
		float a = 1.0f / (alpha * absTanTheta);
		if (a >= 1.6f) return 0.0f;
		return (1.0f - 1.259f * a + 0.396f * a * a) / (3.535f * a + 2.181f * a * a);
	}
	static float Gggx(Vec3 wi, Vec3 wo, float alpha)
	{
		return 1.0f / (1.0f + lambdaGGX(wi, alpha) + lambdaGGX(wo, alpha));
	}
	static float fresnelDielectric(float cosThetaI, float iorExt, float iorInt)
	{
		cosThetaI = fabsf(cosThetaI);
		float etaI = iorExt;
		float etaT = iorInt;
		float sinThetaI = sqrtf(std::max(0.0f, 1.0f - cosThetaI * cosThetaI));
		float eta = etaI / etaT;
		float sinThetaT = eta * sinThetaI;
		if (sinThetaT >= 1.0f) {
			return 1.0f;
		}

		float cosThetaT = sqrtf(std::max(0.0f, 1.0f - sinThetaT * sinThetaT));
		float Rs = (etaI * cosThetaI - etaT * cosThetaT) / (etaI * cosThetaI + etaT * cosThetaT);
		float Rp = (etaT * cosThetaI - etaI * cosThetaT) / (etaT * cosThetaI + etaI * cosThetaT);
		return 0.5f * (Rs * Rs + Rp * Rp);
	}
	static Colour fresnelCondutor(float cosTheta, Colour ior, Colour k)
	{
		cosTheta = clamp(cosTheta, 0.0f, 1.0f);
		Colour cos2 = Colour(cosTheta * cosTheta, cosTheta * cosTheta, cosTheta * cosTheta);
		Colour sin2 = Colour(1.0f, 1.0f, 1.0f) - cos2;
		Colour eta2 = ior * ior;
		Colour k2 = k * k;
		Colour t0 = eta2 - k2 - sin2;
		Colour a2plusb2 = (t0 * t0 + Colour(4.0f, 4.0f, 4.0f) * eta2 * k2).sqrt();
		Colour a = ((a2plusb2 + t0) * 0.5f).sqrt();
		Colour temp = a * Colour(2.0f, 2.0f, 2.0f) * Colour(cosTheta, cosTheta, cosTheta);
		Colour Rs = ((a2plusb2 + cos2) - temp) / ((a2plusb2 + cos2) + temp);
		Colour Rp = Rs * ((cos2 * a2plusb2 + sin2 * sin2) - (temp * sin2)) /((cos2 * a2plusb2 + sin2 * sin2) + (temp * sin2));
		return (Rs + Rp) * 0.5f;
	}
};

class BSDF
{
public:
	Colour emission;
	virtual Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf) = 0;
	virtual Colour evaluate(const ShadingData& shadingData, const Vec3& wi) = 0;
	virtual float PDF(const ShadingData& shadingData, const Vec3& wi) = 0;
	virtual bool isPureSpecular() = 0;
	virtual bool isTwoSided() = 0;
	bool isLight()
	{
		return emission.Lum() > 0 ? true : false;
	}
	void addLight(Colour _emission)
	{
		emission = _emission;
	}
	Colour emit(const ShadingData& shadingData, const Vec3& wi)
	{
		return emission;
	}
	virtual float mask(const ShadingData& shadingData) = 0;
};

class DiffuseBSDF : public BSDF
{
public:
	Texture* albedo;
	DiffuseBSDF() = default;
	DiffuseBSDF(Texture* _albedo)
	{
		albedo = _albedo;
	}
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf)
	{
		Vec3 wi_local = SamplingDistributions::cosineSampleHemisphere(sampler->next(), sampler->next());
		pdf = SamplingDistributions::cosineHemispherePDF(wi_local);
		reflectedColour = albedo->sample(shadingData.tu, shadingData.tv) / M_PI;
		Vec3 wi = shadingData.frame.toWorld(wi_local);
		return wi;
	}
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		return albedo->sample(shadingData.tu, shadingData.tv) / M_PI;
	}
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		Vec3 wi_local = shadingData.frame.toLocal(wi);
		return SamplingDistributions::cosineHemispherePDF(wi_local);
	}
	bool isPureSpecular()
	{
		return false;
	}
	bool isTwoSided()
	{
		return true;
	}
	float mask(const ShadingData& shadingData)
	{
		return albedo->sampleAlpha(shadingData.tu, shadingData.tv);
	}
};

class MirrorBSDF : public BSDF
{
public:
	Texture* albedo;
	MirrorBSDF() = default;
	MirrorBSDF(Texture* _albedo)
	{
		albedo = _albedo;
	}
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf)
	{
		Vec3 wo = shadingData.wo;
		Vec3 n = shadingData.sNormal;
		Vec3 wi = ShadingHelper::reflect(wo, n); 

		pdf = 1.0f;
		reflectedColour = albedo->sample(shadingData.tu, shadingData.tv);
		return wi;
	}
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		return Colour(0.0f, 0.0f, 0.0f);
	}
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		return 0.0f;
	}
	bool isPureSpecular()
	{
		return true;
	}
	bool isTwoSided()
	{
		return true;
	}
	float mask(const ShadingData& shadingData)
	{
		return albedo->sampleAlpha(shadingData.tu, shadingData.tv);
	}
};

class ConductorBSDF : public BSDF
{
public:
	Texture* albedo;
	Colour eta;
	Colour k;
	float alpha;
	ConductorBSDF() = default;
	ConductorBSDF(Texture* _albedo, Colour _eta, Colour _k, float roughness)
	{
		albedo = _albedo;
		eta = _eta;
		k = _k;
		alpha = 1.62142f * sqrtf(roughness);
	}
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf)
	{
		Vec3 wo = shadingData.frame.toLocal(shadingData.wo);
		Vec3 h = SamplingDistributions::sampleGGX(wo, alpha, sampler->next(), sampler->next());
		Vec3 wi = ShadingHelper::reflect(wo, h);

		if (wi.z <= 0.0f) {
			pdf = 0.0f;
			reflectedColour = Colour(0.0f, 0.0f, 0.0f);
			return shadingData.wo; 
		}

		float D = ShadingHelper::Dggx(h, alpha);
		float G = ShadingHelper::Gggx(wi, wo, alpha);
		Colour F = ShadingHelper::fresnelCondutor(Dot(wi, h), eta, k);
		float denom = 4.0f * std::max(1e-6f, fabsf(wo.z) * fabsf(wi.z));
		reflectedColour = albedo->sample(shadingData.tu, shadingData.tv) * F * D * G / denom;
		float Jh = 1.0f / (4.0f * std::max(1e-6f, Dot(wo, h)));
		pdf = D * fabsf(h.z) * Jh;

		return shadingData.frame.toWorld(wi);
	}
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		Vec3 wo = shadingData.frame.toLocal(shadingData.wo);
		Vec3 localWi = shadingData.frame.toLocal(wi);
		if (localWi.z <= 1e-4f) {
			return Colour(0.0f, 0.0f, 0.0f);
		}

		Vec3 h = (wo + localWi).normalize();
		float D = ShadingHelper::Dggx(h, alpha);
		float G = ShadingHelper::Gggx(localWi, wo, alpha);
		Colour F = ShadingHelper::fresnelCondutor(Dot(localWi, h), eta, k);
		float denom = 4.0f * std::max(1e-6f, fabsf(wi.z) * fabsf(wo.z));
		return albedo->sample(shadingData.tu, shadingData.tv) * F * D * G / denom;
	}
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		Vec3 wo = shadingData.frame.toLocal(shadingData.wo);
		Vec3 localWi = shadingData.frame.toLocal(wi);
		if (localWi.z <= 0.0f) return 0.0f;
		Vec3 h = (wo + localWi).normalize();
		float Jh = 1.0f / (4.0f * std::max(1e-6f, Dot(wo, h)));
		return ShadingHelper::Dggx(h, alpha) * fabsf(h.z) * Jh;
	}
	bool isPureSpecular()
	{
		return false;
	}
	bool isTwoSided()
	{
		return true;
	}
	float mask(const ShadingData& shadingData)
	{
		return albedo->sampleAlpha(shadingData.tu, shadingData.tv);
	}
};

class GlassBSDF : public BSDF
{
public:
	Texture* albedo;
	float intIOR;
	float extIOR;
	GlassBSDF() = default;
	GlassBSDF(Texture* _albedo, float _intIOR, float _extIOR)
	{
		albedo = _albedo;
		intIOR = _intIOR;
		extIOR = _extIOR;
	}

	inline Vec3 refract(const Vec3& I, float cosThetaI, float eta)
	{
		float sinThetaI2 = std::max(0.0f, 1.0f - cosThetaI * cosThetaI);
		float sinThetaT2 = eta * eta * sinThetaI2;
		float cosThetaT = sqrtf(std::max(0.0f, 1.0f - sinThetaT2));

		Vec3 T;
		T.x = eta * (-I.x);
		T.y = eta * (-I.y);
		if (I.z >= 0.0f)
			T.z = -cosThetaT;
		else
			T.z = cosThetaT;

		return T;
	}


	Vec3 sample(const ShadingData& shadingData, Sampler* sampler,
		Colour& reflectedColour, float& pdf)
	{
		Colour baseColor = albedo->sample(shadingData.tu, shadingData.tv);

		Vec3 woWorld = shadingData.wo;
		Vec3 woLocal = shadingData.frame.toLocal(woWorld);

		float cosThetaO = woLocal.z;
		if (fabsf(cosThetaO) < 1e-6f) {
			pdf = 0.f;
			reflectedColour = Colour(0.f, 0.f, 0.f);
			return Vec3(0.f, 0.f, 0.f);
		}

		bool entering = (cosThetaO > 0.0f);
		float n1 = entering ? extIOR : intIOR;
		float n2 = entering ? intIOR : extIOR;
		float eta = n1 / n2;

		float F = ShadingHelper::fresnelDielectric(cosThetaO, n1, n2);

		float randVal = sampler->next();

		Vec3 wiLocal(0.f, 0.f, 0.f);
		if (randVal < F) {
			Vec3 I = -woLocal;
			float dotVal = I.z;
			wiLocal = Vec3(I.x - 2.0f * dotVal * 0.f,
				I.y - 2.0f * dotVal * 0.f,
				I.z - 2.0f * dotVal * 1.f);
			pdf = F;

			reflectedColour = baseColor;
		}
		else {
			float cosThetaI = woLocal.z;
			Vec3 I = -woLocal;
			float sinThetaI2 = std::max(0.0f, 1.0f - cosThetaI * cosThetaI);
			float eta2 = eta * eta;
			float sinThetaT2 = eta2 * sinThetaI2;
			if (sinThetaT2 >= 1.0f) {
				float dotVal = I.z;
				wiLocal.z = I.z;
				pdf = 1.f;
				reflectedColour = baseColor;
			}
			else {
				wiLocal = refract(I, cosThetaI, eta);
				pdf = 1.0f - F;
				reflectedColour = baseColor;
			}
		}

		wiLocal.normalize();
		Vec3 wiWorld = shadingData.frame.toWorld(wiLocal);

		return wiWorld;
	}

	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		return Colour(0.0f, 0.0f, 0.0f);
	}
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		return 0.0f;
	}
	bool isPureSpecular()
	{
		return true;
	}
	bool isTwoSided()
	{
		return false;
	}
	float mask(const ShadingData& shadingData)
	{
		return albedo->sampleAlpha(shadingData.tu, shadingData.tv);
	}
};

class DielectricBSDF : public BSDF
{
public:
	Texture* albedo;
	float intIOR;
	float extIOR;
	float alpha;
	DielectricBSDF() = default;
	DielectricBSDF(Texture* _albedo, float _intIOR, float _extIOR, float roughness)
	{
		albedo = _albedo;
		intIOR = _intIOR;
		extIOR = _extIOR;
		alpha = 1.62142f * sqrtf(roughness);
	}
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf)
	{
		Vec3 wo = shadingData.frame.toLocal(shadingData.wo);
		bool entering = wo.z > 0.0f;
		float etaI = entering ? extIOR : intIOR;
		float etaT = entering ? intIOR : extIOR;
		float eta = etaI / etaT;
		Vec3 h = SamplingDistributions::sampleGGX(wo, alpha, sampler->next(), sampler->next());
		float cosTheta = Dot(wo, h);
		float F = ShadingHelper::fresnelDielectric(cosTheta, etaI, etaT);
		bool isReflect = (sampler->next() < F);

		if (isReflect)
		{
			Vec3 wi = ShadingHelper::reflect(wo, h);
			if (wi.z <= 0.0f) {
				pdf = 0.0f;
				return shadingData.wo;
			}

			float D = ShadingHelper::Dggx(h, alpha);
			float G = ShadingHelper::Gggx(wi, wo, alpha);
			pdf = F * D * fabsf(h.z) / (4.0f * fabsf(Dot(wo, h)));

			reflectedColour = albedo->sample(shadingData.tu, shadingData.tv) * F * D * G /
				(4.0f * fabsf(wi.z) * fabsf(wo.z));
			return shadingData.frame.toWorld(wi);
		}
		else
		{
			float cosThetaT2 = 1.0f - eta * eta * (1.0f - cosTheta * cosTheta);
			if (cosThetaT2 <= 0.0f)
			{
				pdf = 0.0f;
				return shadingData.wo;
			}

			float cosThetaT = sqrtf(cosThetaT2);
			Vec3 wt = -wo * eta + h * (eta * cosTheta - cosThetaT);
			if (wt.z <= 0.0f) {
				pdf = 0.0f;
				return shadingData.wo;
			}

			float D = ShadingHelper::Dggx(h, alpha);
			float G = ShadingHelper::Gggx(wt, wo, alpha);
			float sqrtDenom = Dot(wo, h) + eta * Dot(wt, h);
			float factor = (Dot(wt, h) * Dot(wt, h)) / (Dot(wo, h) + eta * Dot(wt, h));
			pdf = (1.0f - F) * D * fabsf(h.z) * eta * eta * fabsf(Dot(wt, h)) /
				(sqrtDenom * sqrtDenom);

			reflectedColour = albedo->sample(shadingData.tu, shadingData.tv) * (1.0f - F) * D * G * eta * eta *
				fabsf(Dot(wt, h)) * fabsf(Dot(wo, h)) /
				(fabsf(wo.z) * fabsf(wt.z) * sqrtDenom * sqrtDenom);

			return shadingData.frame.toWorld(wt);
		}
	}
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		return Colour(0.0f, 0.0f, 0.0f);
	}
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		return 0.0f;
	}
	bool isPureSpecular()
	{
		return false;
	}
	bool isTwoSided()
	{
		return false;
	}
	float mask(const ShadingData& shadingData)
	{
		return albedo->sampleAlpha(shadingData.tu, shadingData.tv);
	}
};

class OrenNayarBSDF : public BSDF
{
public:
	Texture* albedo;
	float sigma;
	OrenNayarBSDF() = default;
	OrenNayarBSDF(Texture* _albedo, float _sigma)
	{
		albedo = _albedo;
		sigma = _sigma;
	}
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf)
	{
		Vec3 wiLocal = SamplingDistributions::cosineSampleHemisphere(
			sampler->next(), sampler->next());
		pdf = SamplingDistributions::cosineHemispherePDF(wiLocal);
		reflectedColour = albedo->sample(shadingData.tu, shadingData.tv);
		reflectedColour = reflectedColour / M_PI;
		Vec3 wi = shadingData.frame.toWorld(wiLocal);
		return wi;
	}
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		Vec3 woLocal = shadingData.frame.toLocal(shadingData.wo);
		Vec3 wiLocal = shadingData.frame.toLocal(wi);

		if (woLocal.z <= 0.0f || wiLocal.z <= 0.0f) return Colour(0.0f, 0.0f, 0.0f);

		float sigma2 = sigma * sigma;
		float A = 1.0f - (0.5f * sigma2 / (sigma2 + 0.33f));
		float B = 0.45f * sigma2 / (sigma2 + 0.09f);
		float cosThetaI = wiLocal.z;
		float cosThetaO = woLocal.z;
		float alpha = std::max(acosf(cosThetaI), acosf(cosThetaO));
		float beta = std::min(acosf(cosThetaI), acosf(cosThetaO));
		float phiI = atan2f(wiLocal.y, wiLocal.x);
		float phiO = atan2f(woLocal.y, woLocal.x);
		float deltaPhi = phiI - phiO;
		float cosDeltaPhi = cosf(deltaPhi);
		float sinAlpha = sinf(alpha);
		float tanBeta = tanf(beta);

		float C = A + B * std::max(0.0f, cosDeltaPhi) * sinAlpha * tanBeta;
		Colour base = albedo->sample(shadingData.tu, shadingData.tv);
		base = base * (C / M_PI);
		return base;
	}
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		Vec3 wiLocal = shadingData.frame.toLocal(wi);
		return SamplingDistributions::cosineHemispherePDF(wiLocal);
	}
	bool isPureSpecular()
	{
		return false;
	}
	bool isTwoSided()
	{
		return true;
	}
	float mask(const ShadingData& shadingData)
	{
		return albedo->sampleAlpha(shadingData.tu, shadingData.tv);
	}
};

class PlasticBSDF : public BSDF
{
public:
	Texture* albedo;
	float intIOR;
	float extIOR;
	float alpha;
	PlasticBSDF() = default;
	PlasticBSDF(Texture* _albedo, float _intIOR, float _extIOR, float roughness)
	{
		albedo = _albedo;
		intIOR = _intIOR;
		extIOR = _extIOR;
		alpha = 1.62142f * sqrtf(roughness);
	}

	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf)
	{
		Vec3 wo = shadingData.frame.toLocal(shadingData.wo);
		float cosThetaO = fabsf(wo.z);
		float F = ShadingHelper::fresnelDielectric(cosThetaO, extIOR, intIOR);

		if (sampler->next() < F)
		{
			Vec3 h = SamplingDistributions::sampleGGX(wo, alpha, sampler->next(), sampler->next());
			Vec3 wi = ShadingHelper::reflect(wo, h);
			if (wi.z <= 0.0f) {
				pdf = 0.0f;
				reflectedColour = Colour(0.0f, 0.0f, 0.0f);
				return shadingData.frame.toWorld(Vec3(0, 0, 1));
			}
			float D = ShadingHelper::Dggx(h, alpha);
			float G = ShadingHelper::Gggx(wi, wo, alpha);
			Colour albedoVal = albedo->sample(shadingData.tu, shadingData.tv);
			reflectedColour = albedoVal * D * G * F / (4.0f * fabsf(wi.z) * fabsf(wo.z));
			pdf = D * fabsf(h.z) / (4.0f * Dot(wo, h));
			return shadingData.frame.toWorld(wi);
		}
		else
		{
			Vec3 wiLocal = SamplingDistributions::cosineSampleHemisphere(sampler->next(), sampler->next());
			reflectedColour = albedo->sample(shadingData.tu, shadingData.tv) * (1.0f - F) / M_PI;
			pdf = (1.0f - F) * SamplingDistributions::cosineHemispherePDF(wiLocal);
			return shadingData.frame.toWorld(wiLocal);
		}
	}
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		Vec3 wo = shadingData.frame.toLocal(shadingData.wo);
		Vec3 wiLocal = shadingData.frame.toLocal(wi);
		if (wiLocal.z <= 0.0f) return Colour();
		float cosThetaO = fabsf(wo.z);
		float F = ShadingHelper::fresnelDielectric(cosThetaO, extIOR, intIOR);
		Colour albedoVal = albedo->sample(shadingData.tu, shadingData.tv);
		Vec3 h = (wo + wiLocal).normalize();
		float D = ShadingHelper::Dggx(h, alpha);
		float G = ShadingHelper::Gggx(wiLocal, wo, alpha);
		Colour spec = albedoVal * D * G * F / (4.0f * fabsf(wiLocal.z) * fabsf(wo.z));
		Colour diff = albedoVal * (1.0f - F) / M_PI;
		return spec + diff;
	}
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		Vec3 wo = shadingData.frame.toLocal(shadingData.wo);
		Vec3 wiLocal = shadingData.frame.toLocal(wi);
		if (wiLocal.z <= 0.0f) return 0.0f;
		float cosThetaO = fabsf(wo.z);
		float F = ShadingHelper::fresnelDielectric(cosThetaO, extIOR, intIOR);
		Vec3 h = (wo + wiLocal).normalize();
		float D = ShadingHelper::Dggx(h, alpha);
		float pdfSpec = D * fabsf(h.z) / (4.0f * Dot(wo, h));
		float pdfDiff = SamplingDistributions::cosineHemispherePDF(wiLocal);
		return F * pdfSpec + (1.0f - F) * pdfDiff;
	}
	bool isPureSpecular()
	{
		return false;
	}
	bool isTwoSided()
	{
		return true;
	}
	float mask(const ShadingData& shadingData)
	{
		return albedo->sampleAlpha(shadingData.tu, shadingData.tv);
	}
};

class LayeredBSDF : public BSDF
{
public:
	BSDF* base;
	Colour sigmaa;
	float thickness;
	float intIOR;
	float extIOR;
	LayeredBSDF() = default;
	LayeredBSDF(BSDF* _base, Colour _sigmaa, float _thickness, float _intIOR, float _extIOR)
	{
		base = _base;
		sigmaa = _sigmaa;
		thickness = _thickness;
		intIOR = _intIOR;
		extIOR = _extIOR;
	}
	Vec3 sample(const ShadingData& shadingData, Sampler* sampler, Colour& reflectedColour, float& pdf)
	{
		// Add code to include layered sampling
		return base->sample(shadingData, sampler, reflectedColour, pdf);
	}
	Colour evaluate(const ShadingData& shadingData, const Vec3& wi)
	{
		// Add code for evaluation of layer
		return base->evaluate(shadingData, wi);
	}
	float PDF(const ShadingData& shadingData, const Vec3& wi)
	{
		// Add code to include PDF for sampling layered BSDF
		return base->PDF(shadingData, wi);
	}
	bool isPureSpecular()
	{
		return base->isPureSpecular();
	}
	bool isTwoSided()
	{
		return true;
	}
	float mask(const ShadingData& shadingData)
	{
		return base->mask(shadingData);
	}
};