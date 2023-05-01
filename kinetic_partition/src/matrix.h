#pragma once
#pragma warning(push, 0)
#pragma warning(disable: 4251)
#include "gdal.h"
#include "gdal_priv.h"
#include "cpl_port.h"
#pragma warning(pop)
#include "defs.h"
#include <algorithm>
#include <list>
#include <utility>

using std::list;
using std::pair;


template<typename T>
class Matrix {
public:
	Matrix()
	{
		rows = 0;
		cols = 0;
		channels = 0;
		
		data = NULL;
	}


	Matrix(Matrix & M)
	{
		rows = M.rows;
		cols = M.cols;
		channels = M.channels;
		data = NULL;

		if (rows > 0 && cols > 0 && channels > 0) {
			data = new T*[rows];
			for (uint i = 0; i < rows; i++) {
				data[i] = new T[channels * cols];
				for (uint j = 0; j < cols; j++) {
					for (uint c = 0; c < channels; c++) {
						data[i][j * channels + c] = M(i, j, c);
					}
				}
			}
		}
	}


	Matrix(const uint _rows, const uint _cols, const uint _channels)
	{
		rows = _rows;
		cols = _cols;
		channels = _channels;
		
		data = new T*[rows];
		for (uint i = 0 ; i < rows ; i++) {
            data[i] = new T[channels * cols]();
		}
	}

	Matrix(const uint _rows, const uint _cols, const uint _channels, const T value)
	{
		rows = _rows;
		cols = _cols;
		channels = _channels;

		data = new T*[rows];
		for (uint i = 0; i < rows; i++) {
			int chan_cols = channels * cols;
			T* row_data = new T[chan_cols]();
			for (uint j = 0; j < chan_cols; ++j) row_data[j] = value;
			data[i] = row_data;
		}
	}


    Matrix(GDALDataset* im_dataset)
	{
		// Gets the size and the depth of the image
		// and allocates a chunk with the appropriate length
		cols = im_dataset->GetRasterXSize();
		rows = im_dataset->GetRasterYSize();
		channels = im_dataset->GetRasterCount();

		data = new T*[rows];
		for (uint i = 0 ; i < rows ; i++) data[i] = new T[cols * channels];

		//GByte** any_raster = new GByte*[rows];
		T** im_values = new T*[rows];
		for (uint c = 0 ; c < channels ; c++) {

			// Sequentially loads the raster bands
			GDALRasterBand* im_band = im_dataset->GetRasterBand(c + 1);
			
			for (uint i = 0; i < rows; i++) {
				//any_raster[i] = (GByte*)CPLMalloc(cols * sizeof(GByte));
				//im_band->RasterIO(GF_Read, 0, i, cols, 1, any_raster[i], cols, 1, GDT_Byte, 0, 0);
				im_values[i] = (T*)CPLMalloc(cols * sizeof(T));
				if (std::is_same<T, uchar>::value) {
					im_band->RasterIO(GF_Read, 0, i, cols, 1, im_values[i], cols, 1, GDT_Byte, 0, 0);
				}
				else {
					im_band->RasterIO(GF_Read, 0,  i, cols, 1, im_values[i], cols, 1, GDT_Float32, 0, 0);
				}
			}

			// Copies data
			for (uint i = 0 ; i < rows ; i++) {
				for (uint j = 0 ; j < cols ; j++) {
					//data[i][channels * j + c] = static_cast<T>(any_raster[i][j]);
					data[i][channels * j + c] = static_cast<T>(im_values[i][j]);
				}
			}
			//for (uint i = 0; i < rows; i++) CPLFree(any_raster[i]);
			for (uint i = 0; i < rows; i++) CPLFree(im_values[i]);
		}

		//delete[] any_raster;
		delete[] im_values;
	}


	Matrix(GDALDataset* im_dataset, uint size_i, uint size_j, uint offset_i, uint offset_j)
	{
		// Gets the size and the depth of the image
		// and allocates a chunk with the appropriate length
		cols = size_j;
		rows = size_i;
		channels = im_dataset->GetRasterCount();

		data = new T*[rows];
		for (uint i = 0; i < rows; i++) data[i] = new T[cols * channels];

		T** im_values = new T*[rows];
		for (uint c = 0; c < channels; c++) {

			// Sequentially loads the raster bands
			GDALRasterBand* im_band = im_dataset->GetRasterBand(c + 1);
			for (uint i = 0; i < rows; i++) {
				im_values[i] = (T*)CPLMalloc(cols * sizeof(T));
				if (std::is_same<T, uchar>::value) {
					im_band->RasterIO(GF_Read, offset_j, offset_i + i, cols, 1, im_values[i], cols, 1, GDT_Byte, 0, 0);
				} else {
					im_band->RasterIO(GF_Read, offset_j, offset_i + i, cols, 1, im_values[i], cols, 1, GDT_Float32, 0, 0);
				}
			}

			// Copies data
			for (uint i = 0; i < rows; i++) {
				for (uint j = 0; j < cols; j++) {
					data[i][channels * j + c] = im_values[i][j];
				}
			}

			for (uint i = 0; i < rows; i++) CPLFree(im_values[i]);
		}
		delete[] im_values;
    }


	Matrix(Matrix<T> & M, uint _rows, uint _cols, uint offset_rows, uint offset_cols)
	{
		// Gets the size and the depth of the image
		// and allocates a chunk with the appropriate length
		cols = _cols;
		rows = _rows;
		channels = M.channels;

		data = new T*[rows];
		for (uint i = 0; i < rows; i++) data[i] = new T[cols * channels];

		for (uint c = 0; c < channels; c++) {
			for (uint i = 0 ; i < rows ; i++) {
				for (uint j = 0 ; j < cols ; j++) {
					data[i][channels * j + c] = M(offset_rows + i, offset_cols + j, c);
				}
			}
		}
	}


	Matrix(Matrix<T> & M, uint _rows, uint _cols, uint offset_rows, uint offset_cols, uint _channels)
	{
		// Gets the size and the depth of the image
		// and allocates a chunk with the appropriate length
		cols = _cols;
		rows = _rows;
		channels = _channels;

		data = new T*[rows];
		for (uint i = 0; i < rows; i++) data[i] = new T[cols * channels];

		for (uint c = 0; c < channels; c++) {
			for (uint i = 0; i < rows; i++) {
				for (uint j = 0; j < cols; j++) {
					data[i][channels * j + c] = M(offset_rows + i, offset_cols + j, c);
				}
			}
		}
	}


    Matrix(Matrix<T> & M, double scale)
    {
        cols = uint(scale * M.cols);
        rows = uint(scale * M.rows);
        channels = M.channels;

        data = new T*[rows];
        for (uint i = 0 ; i < rows ; i++) data[i] = new T[cols * channels];

        if (scale == 1) {
            // Source and destination images have the same size
            // We just copy it
            for (uint c = 0; c < channels; c++) {
                for (uint i = 0; i < rows; i++) {
                    for (uint j = 0; j < cols; j++) {
                        data[i][channels * j + c] = M(i, j, c);
                    }
                }
            }
        } else if (scale < 1) {
            // The source image is larger than the destination image

            Matrix<double> G;
            int radius = int(ceil(1.0 / scale));
            gaussian_matrix(G, radius);

            for (uint c = 0; c < channels; c++) {
                for (uint i = 0; i < rows; i++) {
                    for (uint j = 0; j < cols; j++) {
                        double result = 0;
                        double normalizer = 0;
                        for (int k = -radius ; k <= radius ; k++) {
                            if ((int(int(i) / scale + k) < 0) || (int(int(i) / scale + k) > int(M.rows) - 1)) continue;
                            for (int l = -radius ; l <= radius ; l++) {
                                if ((int(int(j) / scale + l) < 0) || (int(int(j) / scale + l) > int(M.cols) - 1)) continue;
                                result += M(uint(int(i) / scale + k), uint(int(j) / scale + l), c) * G(uint(k + radius), uint(l + radius));
                                normalizer += G(uint(k + radius), uint(l + radius));
                            }
                        }
                        (*this)(i, j, c) = T(result / normalizer);
                    }
                }
            }

        } else if (scale > 1) {
            // The destination image is larger than the source image
            Matrix<double> G;
            int radius = int(ceil(scale));
            gaussian_matrix(G, radius);

            Matrix<T> N (uint(scale * M.rows), uint(scale * M.cols), M.channels);
            for (uint c = 0 ; c < M.channels ; c++) {
                for (uint i = 0 ; i < M.rows ; i++) {
                    uint k_min = uint(scale * i);
                    uint k_max = uint(scale * (i + 1) < N.rows ? scale * (i + 1) : N.rows);
                    for (uint k = k_min ; k < k_max ; k++) {
                        for (uint j = 0 ; j < M.cols ; j++) {
                            uint l_min = uint(scale * j);
                            uint l_max = uint(scale * (j + 1) < N.cols ? scale * (j + 1) : N.cols);
                            for (uint l = l_min ; l <l_max ; l++) {
                                N (k, l, c) = M (i, j, c);
                            }
                        }
                    }
                }
            }

            for (uint c = 0 ; c < channels ; c++) {
                for (uint i = 0 ; i < rows ; i++) {
                    for (uint j = 0 ; j < cols ; j++) {
                        double result = 0;
                        double normalizer = 0;
                        for (int k = -radius ; k <= radius ; k++) {
                            if ((int(i) + k < 0) || (int(i) + k > int(N.rows) - 1)) continue;
                            for (int l = -radius ; l <= radius ; l++) {
                                if ((int(j) + l < 0) || (int(j) + l > int(N.cols) - 1)) continue;
                                result += G(uint(k + radius), uint(l + radius)) * N(i + k, j + l, c);
                                normalizer += G(uint(k + radius), uint(l + radius));
                            }
                        }
                        (*this)(i, j, c) = T(result / normalizer);
                    }
                }
            }
        }
    }


	~Matrix() 
	{
		if (data != NULL) {
			for (uint i = 0; i < rows; i++) {
				delete[] data[i];
			}
			delete[] data;
		}
	}


	T & operator()(uint i, uint j, uint k = 0) 
	{
		return (channels == 1 ? data[i][j] : data[i][channels * j + k]);
	}


	const T & operator()(uint i, uint j, uint k = 0) const
	{
		return (channels == 1 ? data[i][j] : data[i][channels * j + k]);
	}


    Matrix & operator=(const Matrix<T> & M)
    {
        if (this != &M) {
            // If 'this' has a different size from M, we must deallocate and reallocate the data
            if (rows != M.rows || cols != M.cols || channels != M.channels) {
                if (data != NULL) {
                    for (uint i = 0 ; i < rows ; i++) delete[] data[i];
                    delete[] data;
                    data = NULL;
                }
                rows = M.rows;
                cols = M.cols;
                channels = M.channels;
                if (rows > 0 && cols > 0 && channels > 0) {
                    data = new T*[rows];
                    for (uint i = 0 ; i < rows ; i++) {
                        data[i] = new T[channels * cols]();
                    }
                }
            }

            for (uint i = 0 ; i < rows ; i++) {
                memcpy(data[i], M.data[i], channels * cols);
                /*for (uint j = 0 ; j < cols ; j++) {
                    for (uint c = 0 ; c < channels ; c++) {
                        data[i][j * channels + c] = M(i, j, c);
                    }
                }*/
            }
        }
        return *this;
    }


    void set(const T value)
    {
        if (rows == 0 || cols == 0) return;

        for (uint i = 0 ; i < rows ; i++) {
            memset(data[i], value, cols * channels * sizeof(T));
        }
    }


    void release()
    {
        for (uint i = 0 ; i < rows ; i++) {
            delete[] data[i];
        }
        delete[] data;
        data = NULL;

        rows = 0;
        cols = 0;
        channels = 0;
    }


    bool empty()
    {
        return (rows == 0 && cols == 0);
    }

	static void gaussian_matrix_3d(Matrix<double> & G, int radius, double sigma)
	{
		if (!G.empty()) G.release();

		int n = 2 * radius + 1;
		double two_sigma2 = 2 * sigma * sigma;
		double sum = 0;

		G = Matrix<double>(n, n, n);
		for (uint i = 0; i < n; i++) {
			int d_i2 = (i - radius) * (i - radius);
			for (uint j = 0; j < n; j++) {
				int d_j2 = (j - radius) * (j - radius);
				for (uint k = 0; k < n; k++) {
					int d_k2 = (k - radius) * (k - radius);
					G(i, j, k) = exp(-(d_i2 + d_j2 + d_k2) / two_sigma2);
					sum += G(i, j, k);
				}
			}
		}
		for (uint i = 0; i < n; i++) {
			for (uint j = 0; j < n; j++) {
				for (uint k = 0; k < n; k++) {
					G(i, j, k) /= sum;
				}
			}
		}
	}

    static void gaussian_matrix(Matrix<double> & G, int radius, double sigma = 0.3)
    {
        if (!G.empty()) G.release();

        int n = 2 * radius + 1;
        double two_sigma2 = 2 * sigma * sigma;
        double sum = 0;

        G = Matrix<double>(n, n, 1);
        for (uint i = 0 ; i < G.rows ; i++) {
            int d_i2 = (i - radius) * (i - radius);
            for (uint j = 0 ; j < G.cols ; j++) {
                int d_j2 = (j - radius) * (j - radius);
                G(i, j) = exp(-(d_i2 + d_j2) / two_sigma2);
                sum += G(i, j);
            }
        }

        for (uint i = 0 ; i < G.rows ; i++) {
            for (uint j = 0 ; j < G.cols ; j++) {
                G(i, j) /= sum;
            }
        }
    }

    bool write_uchar(std::string & path)
    {
        // We need to load the driver that corresponds to this extension
        int path_size = int(path.size());
        std::string extension;
        for (int i = path_size - 1 ; i >= 0 ; i--) {
            if (path[i] == '.') {
                extension = path.substr(i + 1, path_size - i - 1);
                break;
            }
        }
        GDALDriver* driver = GetGDALDriverManager()->GetDriverByName("GTiff");

        if (driver == NULL) return false;

        // Creates a dataset
        GDALDataset* dataset = driver->Create(path.c_str(), cols, rows, channels, GDT_Byte, NULL);

        // Writes data
        uchar* values = new uchar[cols];
        for (uint c = 0 ; c < channels ; c++) {
            GDALRasterBand* band = dataset->GetRasterBand(c + 1);
            for (uint i = 0 ; i < rows ; i++) {
                for (uint j = 0; j < cols; j++) {
                    values[j] = uchar(jclamp(0, data[i][j * channels + c], 255));
                }
                band->RasterIO(GF_Write, 0, i, cols, 1, values, cols, 1, GDT_Byte, 0, 0);
            }
        }
        delete[] values;

        // Closes dataset
        GDALClose(dataset);
        return true;
    }


    bool write_short(std::string & path)
    {
        // We need to load the driver that corresponds to this extension
        int path_size = int(path.size());
        std::string extension;
        for (int i = path_size - 1; i >= 0; i--) {
            if (path[i] == '.') {
                extension = path.substr(i + 1, path_size - i - 1);
                break;
            }
        }
        GDALDriver* driver = GetGDALDriverManager()->GetDriverByName("GTiff");
        if (driver == NULL) return false;

        // Creates a dataset
        GDALDataset* dataset = driver->Create(path.c_str(), cols, rows, channels, GDT_Int16, NULL);

        // Writes data
        short* values = new short[cols];
        for (uint c = 0; c < channels; c++) {
            GDALRasterBand* band = dataset->GetRasterBand(c + 1);
            for (uint i = 0; i < rows; i++) {
                for (uint j = 0; j < cols; j++) {
                    values[j] = short(data[i][j * channels + c]);
                }
                band->RasterIO(GF_Write, 0, i, cols, 1, values, cols, 1, GDT_Int16, 0, 0);
            }
        }
        delete[] values;

        // Closes dataset
        GDALClose(dataset);
        return true;
    }


    void line(const int y1, const int x1, const int y2, const int x2, const T* values)
	{
		int x = x1, y = y1;
		int dx = x2 - x1;
		int dy = y2 - y1;
		if (dx != 0) {
			if (dx > 0) {
				if (dy != 0) {
					if (dy > 0) {
						if (dx >= dy) {
							int e = dx;
							dx = 2 * e; dy = 2 * dy;
							while (true) {
								for (uint c = 0; c < channels; c++) data[y][channels * x + c] = values[c];
								if (++x == x2) break;
								e -= dy;
								if (e < 0) {
									++y;
									e += dx;
								}
							}
						} else {
							int e = dy;
							dy = 2 * e; dx = 2 * dx;
							while (true) {
								for (uint c = 0; c < channels; c++) data[y][channels * x + c] = values[c];
								if (++y == y2) break;
								e -= dx;
								if (e < 0) {
									++x;
									e += dy;
								}
							}
						}
					} else {
						if (dx >= -dy) {
							int e = dx;
							dx = 2 * e; dy = 2 * dy;
							while (true) {
								for (uint c = 0; c < channels; c++) data[y][channels * x + c] = values[c];
								if (++x == x2) break;
								e += dy;
								if (e < 0) {
									--y;
									e += dx;
								}
							}
						} else {
							int e = dy;
							dy = 2 * e; dx = 2 * dx;
							while (true) {
								for (uint c = 0; c < channels; c++) data[y][channels * x + c] = values[c];
								if (--y == y2) break;
								e += dx;
								if (e > 0) {
									++x;
									e += dy;
								}
							}
						}
					}
					for (uint c = 0; c < channels; c++) data[y2][channels * x2 + c] = values[c];
				} else {
					for (x = x1; x <= x2; x++) {
						for (uint c = 0; c < channels; c++) data[y][channels * x + c] = values[c];
					}
				}
			} else {
				if (dy != 0) {
					if (dy > 0) {
						if (-dx >= dy) {
							int e = dx;
							dx = 2 * e; dy = 2 * dy;
							while (true) {
								for (uint c = 0; c < channels; c++) data[y][channels * x + c] = values[c];
								if (--x == x2) break;
								e += dy;
								if (e >= 0) {
									++y;
									e += dx;
								}
							}
						} else {
							int e = dy;
							dy = 2 * e; dx = 2 * dx;
							while (true) {
								for (uint c = 0; c < channels; c++) data[y][channels * x + c] = values[c];
								if (++y == y2) break;
								e += dx;
								if (e <= 0) {
									--x;
									e += dy;
								}
							}
						}
					} else {
						if (dx <= dy) {
							int e = dx;
							dx = 2 * e; dy = 2 * dy;
							while (true) {
								for (uint c = 0; c < channels; c++) data[y][channels * x + c] = values[c];
								if (--x == x2) break;
								e -= dy;
								if (e >= 0) {
									--y;
									e += dx;
								}
							}
						} else {
							int e = dy;
							dy = 2 * e; dx = 2 * dx;
							while (true) {
								for (uint c = 0; c < channels; c++) data[y][channels * x + c] = values[c];
								if (--y == y2) break;
								e -= dx;
								if (e >= 0) {
									--x;
									e += dy;
								}
							}
						}
					}
					for (uint c = 0; c < channels; c++) data[y2][channels * x2 + c] = values[c];
				} else {
					for (x = x2; x <= x1; x++) {
						for (uint c = 0; c < channels; c++) data[y][channels * x + c] = values[c];
					}
				}
			}
		} else {
			if (dy > 0) {
				for (y = y1; y <= y2; y++) {
					for (uint c = 0; c < channels; c++) data[y][channels * x + c] = values[c];
				}
			} else {
				for (y = y2; y <= y1; y++) {
					for (uint c = 0; c < channels; c++) data[y][channels * x + c] = values[c];
				}
			}
		}
	}


    void line(const int y1, const int x1, const int y2, const int x2, list<pair<int, int> > & pixels)
    {
        pixels.clear();
        int x = x1, y = y1;
        int dx = x2 - x1;
        int dy = y2 - y1;
        if (dx != 0) {
            if (dx > 0) {
                if (dy != 0) {
                    if (dy > 0) {
                        if (dx >= dy) {
                            int e = dx;
                            dx = 2 * e; dy = 2 * dy;
                            while (true) {
                                pixels.push_back(std::make_pair(y, x));
                                if (++x == x2) break;
                                e -= dy;
                                if (e < 0) {
                                    ++y;
                                    e += dx;
                                }
                            }
                        } else {
                            int e = dy;
                            dy = 2 * e; dx = 2 * dx;
                            while (true) {
                                pixels.push_back(std::make_pair(y, x));
                                if (++y == y2) break;
                                e -= dx;
                                if (e < 0) {
                                    ++x;
                                    e += dy;
                                }
                            }
                        }
                    } else {
                        if (dx >= -dy) {
                            int e = dx;
                            dx = 2 * e; dy = 2 * dy;
                            while (true) {
                                pixels.push_back(std::make_pair(y, x));
                                if (++x == x2) break;
                                e += dy;
                                if (e < 0) {
                                    --y;
                                    e += dx;
                                }
                            }
                        } else {
                            int e = dy;
                            dy = 2 * e; dx = 2 * dx;
                            while (true) {
                                pixels.push_back(std::make_pair(y, x));
                                if (--y == y2) break;
                                e += dx;
                                if (e > 0) {
                                    ++x;
                                    e += dy;
                                }
                            }
                        }
                    }
                    pixels.push_back(std::make_pair(y2, x2));
                } else {
                    for (x = x1; x <= x2; x++) {
                        pixels.push_back(std::make_pair(y, x));
                    }
                }
            } else {
                if (dy != 0) {
                    if (dy > 0) {
                        if (-dx >= dy) {
                            int e = dx;
                            dx = 2 * e; dy = 2 * dy;
                            while (true) {
                                pixels.push_back(std::make_pair(y, x));
                                if (--x == x2) break;
                                e += dy;
                                if (e >= 0) {
                                    ++y;
                                    e += dx;
                                }
                            }
                        } else {
                            int e = dy;
                            dy = 2 * e; dx = 2 * dx;
                            while (true) {
                                pixels.push_back(std::make_pair(y, x));
                                if (++y == y2) break;
                                e += dx;
                                if (e <= 0) {
                                    --x;
                                    e += dy;
                                }
                            }
                        }
                    } else {
                        if (dx <= dy) {
                            int e = dx;
                            dx = 2 * e; dy = 2 * dy;
                            while (true) {
                                pixels.push_back(std::make_pair(y, x));
                                if (--x == x2) break;
                                e -= dy;
                                if (e >= 0) {
                                    --y;
                                    e += dx;
                                }
                            }
                        } else {
                            int e = dy;
                            dy = 2 * e; dx = 2 * dx;
                            while (true) {
                                pixels.push_back(std::make_pair(y, x));
                                if (--y == y2) break;
                                e -= dx;
                                if (e >= 0) {
                                    --x;
                                    e += dy;
                                }
                            }
                        }
                    }
                    pixels.push_back(std::make_pair(y2, x2));
                } else {
                    for (x = x2; x <= x1; x++) {
                        pixels.push_back(std::make_pair(y, x));
                    }
                }
            }
        } else {
            if (dy > 0) {
                for (y = y1; y <= y2; y++) {
                    pixels.push_back(std::make_pair(y, x));
                }
            } else {
                for (y = y2; y <= y1; y++) {
                    pixels.push_back(std::make_pair(y, x));
                }
            }
        }
    }

protected:
	T** data;

public:
	uint rows;
	uint cols;
	uint channels;
};
