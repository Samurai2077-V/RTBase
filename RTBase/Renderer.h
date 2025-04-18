#pragma once

#include "Core.h"
#include "Sampling.h"
#include "Geometry.h"
#include "Imaging.h"
#include "Materials.h"
#include "Lights.h"
#include "Scene.h"
#include "GamesEngineeringBase.h"
#include <thread>
#include <functional>
#include <mutex>

struct VPL
{
	ShadingData sData;
	Colour Le;
	bool isLight;

	VPL() {}
	VPL(ShadingData _sData, Colour _Le, bool _isLight)
		: sData(_sData), Le(_Le), isLight(_isLight) {}
};

extern bool isLightTrace = true;
extern bool isPathTrace = true;

class RayTracer
{
public:
	Scene* scene;
	GamesEngineeringBase::Window* canvas;
	Film* film;
	MTRandom* samplers;
	std::thread** threads;
	int numProcs;
	std::vector<VPL> vpls;

	struct Tile {
		int x0, y0, x1, y1;
		int spp; 
		bool converged;
		Colour mean;
		Colour M2; 
		Tile(int _x0, int _y0, int _x1, int _y1)
			: x0(_x0), y0(_y0), x1(_x1), y1(_y1), spp(0), converged(false), mean(0, 0, 0), M2(0, 0, 0) {}
	};

	void init(Scene* _scene, GamesEngineeringBase::Window* _canvas)
	{
		scene = _scene;
		canvas = _canvas;
		film = new Film();
		film->init((unsigned int)scene->camera.width, (unsigned int)scene->camera.height, new BoxFilter());
		SYSTEM_INFO sysInfo;
		GetSystemInfo(&sysInfo);
		numProcs = sysInfo.dwNumberOfProcessors;
		threads = new std::thread * [numProcs];
		samplers = new MTRandom[numProcs];
		clear();
	}
	void clear()
	{
		film->clear();
	}
	Colour computeDirect(const ShadingData& shadingData, Sampler* sampler)
	{
		if (shadingData.bsdf->isPureSpecular())
			return Colour(0.0f, 0.0f, 0.0f);
		Colour result(0.0f, 0.0f, 0.0f);
		float pmf;
		Light* light = scene->sampleLight(sampler, pmf);
		if (light)
		{
			float lightPdf;
			Colour Le;
			Vec3 wi = light->sample(shadingData, sampler, Le, lightPdf);

			if (lightPdf > 1e-6f && pmf > 1e-6f && Le.Lum() > 0.0f)
			{
				bool visible = false;
				if (light->isArea())
				{
					Vec3 toLight = wi - shadingData.x; 
					float dist2 = toLight.lengthSq();
					Vec3 dir = toLight.normalize();
					if (scene->visible(shadingData.x, wi))
					{
						float cos_x = max(Dot(dir, shadingData.sNormal), 0.0f);
						float cos_l = max(-Dot(dir, light->normal(dir)), 0.0f);
						float G = (cos_x * cos_l) / dist2;
						Colour f = shadingData.bsdf->evaluate(shadingData, dir);
						float bsdfPdf = shadingData.bsdf->PDF(shadingData, dir);
						// MIS 
						float w_light = (lightPdf * pmf) * (lightPdf * pmf);
						float w_bsdf = bsdfPdf * bsdfPdf;
						float weight = w_light / (w_light + w_bsdf + 1e-6f);

						result = result + f * Le * G * weight / (lightPdf * pmf);
					}
				}
				else
				{
					// environmentlight
					if (scene->visible(shadingData.x, shadingData.x + wi * 1e4f))
					{
						Colour f = shadingData.bsdf->evaluate(shadingData, wi);
						float cosTerm = max(Dot(wi, shadingData.sNormal), 0.0f);
						float bsdfPdf = shadingData.bsdf->PDF(shadingData, wi);
						// MIS
						float w_light = (lightPdf * pmf) * (lightPdf * pmf);
						float w_bsdf = bsdfPdf * bsdfPdf;
						float weight = w_light / (w_light + w_bsdf + 1e-6f);

						result = result + f * Le * cosTerm * weight / (lightPdf * pmf);
					}
				}
			}
		}

		float bsdfPdf;
		Colour bsdfVal;
		Vec3 wi_bsdf = shadingData.bsdf->sample(shadingData, sampler, bsdfVal, bsdfPdf);

		if (bsdfPdf > 1e-6f && bsdfVal.Lum() > 0.0f)
		{
			float cosTerm = max(Dot(wi_bsdf, shadingData.sNormal), 0.0f);
			if (cosTerm > 0.0f)
			{
				// shadow 
				if (scene->visible(shadingData.x, shadingData.x + wi_bsdf * 1e4f))
				{
					Colour Le(0.0f, 0.0f, 0.0f);
					float lightPdf = 0.0f;
					float pmfLocal = 1.0f; 
					if (scene->background)
					{
						Le = scene->background->evaluate(wi_bsdf);
						lightPdf = scene->background->PDF(shadingData, wi_bsdf);
					}

					if (lightPdf > 1e-6f && Le.Lum() > 0.0f)
					{
						// MIS
						float w_bsdf = bsdfPdf * bsdfPdf;
						float w_light = (lightPdf * pmfLocal) * (lightPdf * pmfLocal);
						float weight = w_bsdf / (w_bsdf + w_light + 1e-6f);

						result = result + bsdfVal * Le * cosTerm * weight / bsdfPdf;
					}
				}
			}
		}
		return result;
	}

	Colour computeDirectInstantRadiosity(ShadingData shadingData)
	{
		if (shadingData.bsdf->isPureSpecular() == true)
		{
			return Colour(0.0f, 0.0f, 0.0f);
		}

		Colour result(0, 0, 0);
		for (auto vpl : vpls) {
			Colour temp(0, 0, 0);
			Vec3 wi(0, 0, 0);
			float GeomtryTermHalfArea = 0;
			wi = vpl.sData.x - shadingData.x;
			float r2Inv = 1.0f / wi.lengthSq();
			wi = wi.normalize();
			float costheta = max(0, wi.dot(shadingData.sNormal));
			float costhetaL = max(0, -wi.dot(vpl.sData.sNormal));
			GeomtryTermHalfArea = costheta * costhetaL * r2Inv;
			if (GeomtryTermHalfArea > 0) {
				if (!scene->visible(shadingData.x, vpl.sData.x)) {
					GeomtryTermHalfArea = 0;
				}
			}
			else {
				GeomtryTermHalfArea = 0;
			}

			if (vpl.isLight) {
				temp = vpl.Le * shadingData.bsdf->evaluate(shadingData, wi) * GeomtryTermHalfArea;
			}
			else {
				temp = vpl.sData.bsdf->evaluate(vpl.sData, wi) * shadingData.bsdf->evaluate(shadingData, wi) * GeomtryTermHalfArea;
			}

			result = result + temp;
		}
		return result;
	}

	const int MAX_DEPTH = 16;

	Colour pathTrace(Ray& r, Colour& pathThroughput, int depth, Sampler* sampler, bool canHitLight = true)
	{
		IntersectionData intersection = scene->bvh->traverse(r, scene->triangles);
		ShadingData shadingData = scene->calculateShadingData(intersection, r);
		if (shadingData.t < FLT_MAX)
		{
			if (shadingData.bsdf->isLight())
			{
				if (canHitLight == true)
				{
					return pathThroughput * shadingData.bsdf->emit(shadingData, shadingData.wo);
				}
				else
				{
					return Colour(0.0f, 0.0f, 0.0f);
				}
			}
			Colour direct = pathThroughput * computeDirect(shadingData, sampler);
			if (depth > MAX_DEPTH)
			{
				return direct;
			}
			float russianRouletteProbability = min(pathThroughput.Lum(), 0.9f);
			if (sampler->next() < russianRouletteProbability)
			{
				pathThroughput = pathThroughput / russianRouletteProbability;
			}
			else
			{
				return direct;
			}

			Colour indirect;
			float pdf;
			Vec3 wi = shadingData.bsdf->sample(shadingData, sampler, indirect, pdf);
			if (pdf < 1e-6f || !std::isfinite(pdf) || indirect.hasNaN())
			{
				return direct;
			}
			float cosTerm = fmaxf(0.0f, Dot(wi, shadingData.sNormal));
			if (!std::isfinite(cosTerm) || !std::isfinite(indirect.r) || !std::isfinite(indirect.g) || !std::isfinite(indirect.b))
			{
				return direct;
			}

			pathThroughput = pathThroughput * indirect * (cosTerm / pdf);
			r.init(shadingData.x + (wi * EPSILON), wi);
			return (direct + pathTrace(r, pathThroughput, depth + 1, sampler, shadingData.bsdf->isPureSpecular()));
		}
		return scene->background->evaluate(r.dir);
	}

	Colour pathTraceHybrid(Ray& r, Colour& pathThroughput, int depth, Sampler* sampler, bool canHitLight = true)
	{
		IntersectionData intersection = scene->bvh->traverse(r, scene->triangles);
		ShadingData shadingData = scene->calculateShadingData(intersection, r);

		if (shadingData.t < FLT_MAX)
		{
			if (shadingData.bsdf->isLight())
			{
				if (canHitLight == true)
				{
					return pathThroughput * shadingData.bsdf->emit(shadingData, shadingData.wo);
				}
				else
				{
					return Colour(0.0f, 0.0f, 0.0f);
				}
			}
			Colour direct = pathThroughput * computeDirect(shadingData, sampler);
			Colour indirectIR = pathThroughput * computeDirectInstantRadiosity(shadingData);
			Colour combinedDirect = direct + indirectIR;

			if (depth > MAX_DEPTH)
			{
				return combinedDirect;
			}

			float rrProb = min(pathThroughput.Lum(), 0.9f);
			if (sampler->next() < rrProb)
			{
				pathThroughput = pathThroughput / rrProb;
			}
			else
			{
				return combinedDirect;
			}
			float pdf;
			Colour bsdfVal;
			Vec3 wi = shadingData.bsdf->sample(shadingData, sampler, bsdfVal, pdf);
			if (pdf < 1e-6f || !std::isfinite(pdf) || bsdfVal.hasNaN())
			{
				return combinedDirect;
			}

			float cosTheta = fmaxf(0.0f, Dot(wi, shadingData.sNormal));
			if (!std::isfinite(cosTheta) || bsdfVal.hasNaN())
			{
				return combinedDirect;
			}
			pathThroughput = pathThroughput * bsdfVal * (cosTheta / pdf);
			r.init(shadingData.x + (wi * EPSILON), wi);
			Colour indirectPath = pathTraceHybrid(r, pathThroughput, depth + 1, sampler, shadingData.bsdf->isPureSpecular());
			return combinedDirect + indirectPath;
		}

		return scene->background->evaluate(r.dir);
	}

	void connectToCamera(Vec3 p, Vec3 n, Colour col) {
		float x, y;
		if (!scene->camera.projectOntoCamera(p, x, y)) return;

		Vec3 pToCam = scene->camera.origin - p;
		if (pToCam.length() <= 0.0f) return;
		float r2Inv = 1.0f / pToCam.lengthSq();
		pToCam = pToCam.normalize();

		float costheta = max(0, -pToCam.dot(scene->camera.viewDirection));
		float costhetaL = max(0, pToCam.dot(n));
		float GeomtryTermHalfArea = costheta * costhetaL * r2Inv;

		if (GeomtryTermHalfArea > 0) {
			if (!scene->visible(p, scene->camera.origin)) {
				return;
			}
			else {
				float cos2 = costheta * costheta;
				float cos4 = cos2 * cos2;
				float we = 1.0f / (scene->camera.Afilm * cos4);
				int idx = (int)y * film->width + (int)x;
				if (!isPathTrace && idx <= film->width * film->height) {
					film->albedoBuffer[idx] = col;
					film->normalBuffer[idx] = Colour(fabsf(n.x), fabsf(n.y), fabsf(n.z));
				}
				film->splat(x, y, col * GeomtryTermHalfArea * we);
			}
		}
		else {
			return;
		}
	}

	void lightTrace(Sampler* sampler) {
		float pmf;
		Light* light = scene->sampleLight(sampler, pmf);
		if (!light) return;

		float pdfPosition = 0.0f;
		Vec3 lightPos = light->samplePositionFromLight(sampler, pdfPosition);
		float pdfDirection = 0.0f;
		Vec3 wi = light->sampleDirectionFromLight(sampler, pdfDirection);
		wi = wi.normalize();
		float pdfTotal = pmf * pdfPosition;
		if (pdfTotal <= 0.0f) return;
		Colour Le = light->evaluate(-wi);
		Colour col = Le / pdfTotal;
		connectToCamera(lightPos, light->normal(-wi), col);
		Ray r(lightPos + (wi * 0.001f), wi);
		Le = Le * wi.dot(light->normal(-wi));
		pdfTotal *= pdfDirection;
		lightTracePath(r, Colour(1.0f, 1.0f, 1.0f), Le / pdfTotal, sampler, 0);
	}

	void lightTracePath(Ray& r, Colour pathThroughput, Colour Le, Sampler* sampler, int depth) {
		IntersectionData intersection = scene->bvh->traverse(r, scene->triangles);
		ShadingData shadingData = scene->calculateShadingData(intersection, r);

		if (shadingData.t < FLT_MAX)
		{
			Vec3 wi1 = scene->camera.origin - shadingData.x;
			wi1 = wi1.normalize();

			connectToCamera(shadingData.x, shadingData.sNormal, (pathThroughput * shadingData.bsdf->evaluate(shadingData, wi1) * Le));

			if (depth > MAX_DEPTH)
			{
				return;
			}
			float russianRouletteProbability = min(pathThroughput.Lum(), 0.9f);
			if (sampler->next() < russianRouletteProbability)
			{
				pathThroughput = pathThroughput / russianRouletteProbability;
			}
			else
			{
				return;
			}

			float pdf;
			Colour bsdf;
			Vec3 wi = shadingData.bsdf->sample(shadingData, sampler, bsdf, pdf);
			wi = wi.normalize();
			pathThroughput = pathThroughput * bsdf * fabsf(Dot(wi, shadingData.sNormal)) / pdf;

			r.init(shadingData.x + (wi * 0.001f), wi);
			lightTracePath(r, pathThroughput, Le, sampler, depth + 1);
		}
	}

	Colour direct(Ray& r, Sampler* sampler)
	{
		IntersectionData intersection = scene->bvh->traverse(r, scene->triangles);
		ShadingData shadingData = scene->calculateShadingData(intersection, r);
		if (shadingData.t < FLT_MAX)
		{
			if (shadingData.bsdf->isLight())
			{
				return shadingData.bsdf->emit(shadingData, shadingData.wo);
			}
			return computeDirect(shadingData, sampler);
		}
		return Colour(0.0f, 0.0f, 0.0f);
	}

	Colour directInstantRadiosity(Ray& r)
	{
		IntersectionData intersection = scene->bvh->traverse(r, scene->triangles);
		ShadingData shadingData = scene->calculateShadingData(intersection, r);
		if (shadingData.t < FLT_MAX)
		{
			if (shadingData.bsdf->isLight())
			{
				return shadingData.bsdf->emit(shadingData, shadingData.wo);
			}
			return computeDirectInstantRadiosity(shadingData);
		}
		return Colour(0.0f, 0.0f, 0.0f);
	}

	Colour albedo(Ray& r)
	{
		IntersectionData intersection = scene->bvh->traverse(r, scene->triangles);
		ShadingData shadingData = scene->calculateShadingData(intersection, r);
		if (shadingData.t < FLT_MAX)
		{
			if (shadingData.bsdf->isLight())
			{
				return shadingData.bsdf->emit(shadingData, shadingData.wo);
			}
			return shadingData.bsdf->evaluate(shadingData, Vec3(0, 1, 0));
		}
		return scene->background->evaluate(r.dir);
	}
	Colour viewNormals(Ray& r)
	{
		IntersectionData intersection = scene->bvh->traverse(r, scene->triangles);
		if (intersection.t < FLT_MAX)
		{
			ShadingData shadingData = scene->calculateShadingData(intersection, r);
			return Colour(fabsf(shadingData.sNormal.x), fabsf(shadingData.sNormal.y), fabsf(shadingData.sNormal.z));
		}
		return Colour(0.0f, 0.0f, 0.0f);
	}

	void renderTileBased(int tileX, int tileY, int tileSize, int filmWidth, int filmHeight) {
		int startX = tileX * tileSize;
		int startY = tileY * tileSize;
		int endX = min(startX + tileSize, filmWidth);
		int endY = min(startY + tileSize, filmHeight);

		for (int y = startY; y < endY; y++) {
			for (int x = startX; x < endX; x++) {
				unsigned char r, g, b;
				film->tonemap(x, y, r, g, b);
				canvas->draw(x, y, r, g, b);
			}
		}
	}

	void STrender() {
		for (unsigned int y = 0; y < film->height; y++)
		{
			for (unsigned int x = 0; x < film->width; x++)
			{
				float px = x + 0.5f;
				float py = y + 0.5f;
				Ray ray = scene->camera.generateRay(px, py);
				Colour throughput(1.0f, 1.0f, 1.0f);  
				Colour col = pathTrace(ray, throughput, 0, samplers);
				film->splat(px, py, col);
				unsigned char r = (unsigned char)(col.r * 255);
				unsigned char g = (unsigned char)(col.g * 255);
				unsigned char b = (unsigned char)(col.b * 255);
				film->tonemap(x, y, r, g, b);
				canvas->draw(x, y, r, g, b);
			}
		}
	}
	void PathTraceMT() {
		const int width = (int)film->width;
		const int height = (int)film->height;
		const int maxSPP = 64;
		const float threshold = 0.001f;
		const int tileSize = 32;
		const int numTilesX = (width + tileSize - 1) / tileSize;
		const int numTilesY = (height + tileSize - 1) / tileSize;
		std::vector<Tile> tiles;

		for (int ty = 0; ty < numTilesY; ty++) {
			for (int tx = 0; tx < numTilesX; tx++) {
				int x0 = tx * tileSize;
				int y0 = ty * tileSize;
				int x1 = min(x0 + tileSize, width);
				int y1 = min(y0 + tileSize, height);
				tiles.emplace_back(x0, y0, x1, y1);
			}
		}

		std::atomic<int> nextTileIndex(0);

		auto worker = [&](int threadID) {
			Sampler* sampler = &samplers[threadID];
			while (true) {
				int tileIndex = nextTileIndex.fetch_add(1);
				if (tileIndex >= tiles.size()) break;
				Tile& tile = tiles[tileIndex];
				if (tile.converged) return; 
				for (int y = tile.y0; y < tile.y1; y++) {
					for (int x = tile.x0; x < tile.x1; x++) {
						float px = x + 0.5f;
						float py = y + 0.5f;
						Ray ray = scene->camera.generateRay(px, py);
						if (tile.spp == 0 && isPathTrace) {
							int idx = y * width + x;
							film->albedoBuffer[idx] = albedo(ray);
							film->normalBuffer[idx] = viewNormals(ray);
						}
						Colour throughput(1.0f, 1.0f, 1.0f);  
						Colour col = pathTraceHybrid(ray, throughput, 0, sampler);
						film->splat(px, py, col);
						tile.spp++;
						Colour delta = col - tile.mean;
						tile.mean = tile.mean + delta * (1.0f / tile.spp);
						tile.M2 = tile.M2 + delta * (col - tile.mean);

						unsigned char r = (unsigned char)(col.r * 255);
						unsigned char g = (unsigned char)(col.g * 255);
						unsigned char b = (unsigned char)(col.b * 255);
						film->tonemap(x, y, r, g, b);
						canvas->draw(x, y, r, g, b);
					}
				}

				if (tile.spp >= 4) {
					Colour var = tile.M2 / (float)(tile.spp - 1);
					float lumVar = var.Lum();
					if (lumVar < threshold || tile.spp >= maxSPP) {
						tile.converged = true;
					}
				}
			}};
		std::vector<std::thread> threadPool;
		for (int i = 0; i < numProcs; i++)
		{
			threadPool.emplace_back(worker, i);
		}

		for (auto& t : threadPool)
		{
			t.join();
		}
	}

	void LightTraceMT() {
		const int filmWidth = film->width;
		const int filmHeight = film->height;
		const int tileSize = 32;
		int numTilesX = (filmWidth + tileSize - 1) / tileSize;
		int numTilesY = (filmHeight + tileSize - 1) / tileSize;
		int totalTiles = numTilesX * numTilesY;
		std::atomic<int> nextTile(0);
		int paths = (filmWidth * filmHeight) / totalTiles;

		std::vector<std::thread> threads;
		threads.reserve(numProcs);
		auto Splat = [&]() {
			while (true) {
				int tileIndex = nextTile.fetch_add(1);
				if (tileIndex >= totalTiles)
					break;
				int tileX = tileIndex % numTilesX;
				int tileY = tileIndex / numTilesX;
				for (int i = 0; i < paths; ++i) {
					lightTrace(samplers);
				}
			}};
		for (int i = 0; i < numProcs; i++) {
			threads.emplace_back(Splat);
		}
		for (auto& t : threads) {
			if (t.joinable())
				t.join();
		}
		threads.clear();
		nextTile = 0; 
		auto Render = [&]() {
			while (true) {
				int tileIndex = nextTile.fetch_add(1);
				if (tileIndex >= totalTiles)
					break;
				int tileX = tileIndex % numTilesX;
				int tileY = tileIndex / numTilesX;
				renderTileBased(tileX, tileY, tileSize, filmWidth, filmHeight);
			}
			};

		for (int i = 0; i < numProcs; i++) {
			threads.emplace_back(Render);
		}
		for (auto& t : threads) {
			if (t.joinable())
				t.join();
		}
	}

	void render()
	{
		traceVPLs(samplers, 5);
		film->incrementSPP();
		if (isLightTrace) LightTraceMT();
		if (isPathTrace) PathTraceMT();
		float* denoised = new float[film->width * film->height * 3];
		Colour* hdrpixels = new Colour[film->width * film->height];
		for (unsigned int i = 0; i < (film->width * film->height); i++)
		{
			hdrpixels[i] = film->film[i] / (float)film->SPP;
		}
		savePNG("GI.png");
		saveHDR("GI.hdr");
		stbi_write_hdr("denoised.hdr", film->width, film->height, 3, denoised);
	}
	int getSPP()
	{
		return film->SPP;
	}
	void saveHDR(std::string filename)
	{
		film->save(filename);
	}
	void savePNG(std::string filename)
	{
		stbi_write_png(filename.c_str(), canvas->getWidth(), canvas->getHeight(), 3, canvas->getBackBuffer(), canvas->getWidth() * 3);
	}
	void traceVPLs(Sampler* sampler, int N_VPLs) {
		vpls.clear();
		vpls.reserve(N_VPLs);
		for (int i = 0; i < N_VPLs; ++i) {
			float pmf;
			Light* light = scene->sampleLight(sampler, pmf);
			if (!light) continue;

			float pdfPosition = 0.0f;
			Vec3 lightPos = light->samplePositionFromLight(sampler, pdfPosition);
			float pdfDirection = 0.0f;
			Vec3 wi = light->sampleDirectionFromLight(sampler, pdfDirection);
			wi = wi.normalize();
			float pdfTotal = pmf * pdfPosition * (float)N_VPLs;

			if (pdfTotal <= 0.0f) continue;

			Colour Le = light->evaluate(-wi);
			Colour col = Le / pdfTotal;
			vpls.emplace_back(VPL(ShadingData(lightPos, light->normal(wi)), col, true));
			Ray r(lightPos + (wi * 0.001f), wi);
			Le = Le * wi.dot(light->normal(-wi));
			pdfTotal *= pdfDirection;
			VPLTracePath(r, Colour(1, 1, 1), Le / pdfTotal, sampler, 0);
		}
	}

	void VPLTracePath(Ray& r, Colour pathThroughput, Colour Le, Sampler* sampler, int depth) {
		IntersectionData intersection = scene->bvh->traverse(r, scene->triangles);
		ShadingData shadingData = scene->calculateShadingData(intersection, r);

		if (shadingData.t < FLT_MAX)
		{
			if (shadingData.bsdf->isLight())
				return;
			else
				vpls.emplace_back(VPL(shadingData, pathThroughput * Le, shadingData.bsdf->isLight()));
			if (depth > MAX_DEPTH) return;
			float russianRouletteProbability = min(pathThroughput.Lum(), 0.9f);
			if (sampler->next() < russianRouletteProbability)
				pathThroughput = pathThroughput / russianRouletteProbability;
			else
				return;

			float pdf;
			Colour bsdf;
			Vec3 wi = shadingData.bsdf->sample(shadingData, sampler, bsdf, pdf);
			wi = wi.normalize();
			pathThroughput = pathThroughput * bsdf * fabsf(Dot(wi, shadingData.sNormal)) / pdf;
			r.init(shadingData.x + (wi * 0.001f), wi);
			VPLTracePath(r, pathThroughput, Le, sampler, depth + 1);
		}
	}
};