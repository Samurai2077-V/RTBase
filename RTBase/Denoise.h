#pragma once
#include <OpenImageDenoise/oidn.hpp>
#include <iostream>

void denoiseWithOIDN(float* color, float* albedo, float* normal, float* output, int width, int height)
{
	try {
		oidn::DeviceRef device = oidn::newDevice();
		device.commit();
		const size_t size = width * height * 3 * sizeof(float); 
		oidn::BufferRef colorBuf = device.newBuffer(size);
		oidn::BufferRef albedoBuf = albedo ? device.newBuffer(size) : nullptr;
		oidn::BufferRef normalBuf = normal ? device.newBuffer(size) : nullptr;
		oidn::BufferRef outputBuf = device.newBuffer(size);

		std::memcpy(colorBuf.getData(), color, size);
		if (albedo) std::memcpy(albedoBuf.getData(), albedo, size);
		if (normal) std::memcpy(normalBuf.getData(), normal, size);

		oidn::FilterRef filter = device.newFilter("RT");
		filter.setImage("color", colorBuf, oidn::Format::Float3, width, height);
		if (albedo) filter.setImage("albedo", albedoBuf, oidn::Format::Float3, width, height);
		if (normal) filter.setImage("normal", normalBuf, oidn::Format::Float3, width, height);
		filter.setImage("output", outputBuf, oidn::Format::Float3, width, height);
		filter.set("hdr", true);
		filter.commit();
		filter.execute();
		std::memcpy(output, outputBuf.getData(), size);
		const char* message;
		if (device.getError(message) != oidn::Error::None) {
			std::cerr << "[OIDN Error] " << message << std::endl;
		}
	}
	catch (const std::exception& e) {
		std::cerr << "OIDN Exception: " << e.what() << std::endl;
	}
}

void saveHDRTonemapped(const std::string& filename, const float* buffer, int width, int height, int spp, float exposure = 1.0f)
{
	float* tonemapped = new float[width * height * 3];
	for (int i = 0; i < width * height; ++i)
	{
		int idx = i * 3;
		float r = buffer[idx + 0] * exposure / spp;
		float g = buffer[idx + 1] * exposure / spp;
		float b = buffer[idx + 2] * exposure / spp;
		r = max(0.0f, r);
		g = max(0.0f, g);
		b = max(0.0f, b);
		tonemapped[idx + 0] = powf(r, 1.0f / 2.2f);
		tonemapped[idx + 1] = powf(g, 1.0f / 2.2f);
		tonemapped[idx + 2] = powf(b, 1.0f / 2.2f);
	}
	stbi_write_hdr(filename.c_str(), width, height, 3, tonemapped);
	delete[] tonemapped;
}