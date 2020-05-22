/* @author : Oğuzhan Yıldırım 240201056*/
// ------------------------------
// Written by Mustafa Ozuysal
// Contact <mustafaozuysal@iyte.edu.tr> for comments and bug reports
// ------------------------------
// Copyright (c) 2018, Mustafa Ozuysal
// All rights reserved.
// ------------------------------
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the copyright holders nor the
//       names of his/its contributors may be used to endorse or promote products
//       derived from this software without specific prior written permission.
// ------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// ------------------------------
#include "image.hpp"

#include <iostream>
#include <vector> 
#include <fstream>
#include <cmath>
#include <cstring>
#include <ostream>
#include <memory>
#include <algorithm>

using std::cerr;
using std::clog;
using std::cos;
using std::endl;
using std::exp;
using std::ifstream;
using std::ios;
using std::memset;
using std::ofstream;
using std::sin;
using std::string;
using std::unique_ptr;
using std::vector;


namespace ceng391 {

Image::Image(int width, int height, int n_channels, int step)
{
        m_width = width;
        m_height = height;
        m_n_channels = n_channels;

        m_step = m_width*m_n_channels;
        if (m_step < step)
                m_step = step;
        m_data = new uchar[m_step*height];
}

Image::~Image()
{
        delete [] m_data;
}

Image* Image::new_gray(int width, int height)
{
        return new Image(width, height, 1);
}

Image* Image::new_rgb(int width, int height)
{
        return new Image(width, height, 3);
}

void Image::set_rect(int x, int y, int width, int height, uchar red, uchar green, uchar blue)
{
        if (x < 0) {
                width += x;
                x = 0;
        }

        if (y < 0) {
                height += y;
                y = 0;
        }

        if (m_n_channels == 1) {
                int value = 0.3*red + 0.59*green + 0.11*blue;
                if (value > 255)
                        value = 255;
                for (int j = y; j < y+height; ++j) {
                        if (j >= m_height)
                                break;
                        uchar* row_data = data(j);
                        for (int i = x; i < x+width; ++i) {
                                if (i >= m_width)
                                        break;
                                row_data[i] = value;
                        }
                }
        } else if (m_n_channels == 3) {
                for (int j = y; j < y+height; ++j) {
                        if (j >= m_height)
                                break;
                        uchar* row_data = data(j);
                        for (int i = x; i < x+width; ++i) {
                                if (i >= m_width)
                                        break;
                                row_data[i*3]     = red;
                                row_data[i*3 + 1] = green;
                                row_data[i*3 + 2] = blue;
                        }
                }
        }
}

void Image::set_rect(int x, int y, int width, int height, uchar value)
{
        if (x < 0) {
                width += x;
                x = 0;
        }

        if (y < 0) {
                height += y;
                y = 0;
        }

        for (int j = y; j < y+height; ++j) {
                if (j >= m_height)
                        break;
                uchar* row_data = data(j);
                for (int i = x; i < x+width; ++i) {
                        if (i >= m_width)
                                break;
                        for (int c = 0; c < m_n_channels; ++c)
                                row_data[i*m_n_channels + c] = value;
                }
        }
}

void Image::to_grayscale()
{
        if (m_n_channels == 1) {
                return;
        } else if (m_n_channels == 3) {
                int new_step = m_width;
                uchar *new_data = new uchar[new_step * m_height];
                for (int y = 0; y < m_height; ++y) {
                        uchar *row_old = m_data + m_step * y;
                        uchar *row_new = new_data + new_step * y;
                        for (int x = 0; x < m_width; ++x) {
                                uchar red = row_old[3*x];
                                uchar green = row_old[3*x + 1];
                                uchar blue = row_old[3*x + 2];
                                int value = 0.3*red + 0.59*green + 0.11*blue;
                                if (value > 255)
                                        value = 255;
                                row_new[x] = value;
                        }
                }

                delete [] m_data;
                m_data = new_data;
                m_step = new_step;
                m_n_channels = 1;
        }
}

void Image::to_rgb()
{
        if (m_n_channels == 3) {
                return;
        } else if (m_n_channels == 1) {
                int new_step = m_width * 3;
                uchar *new_data = new uchar[new_step * m_height];
                for (int y = 0; y < m_height; ++y) {
                        uchar *row_old = m_data + m_step * y;
                        uchar *row_new = new_data + new_step * y;
                        for (int x = 0; x < m_width; ++x) {
                                uchar value = row_old[x];
                                row_new[3*x]     = value;
                                row_new[3*x + 1] = value;
                                row_new[3*x + 2] = value;
                        }
                }

                delete [] m_data;
                m_data = new_data;
                m_step = new_step;
                m_n_channels = 3;
        }
}

bool Image::write_pnm(const std::string& filename) const
{
        ofstream fout;

        string magic_head = "P5";
        string extended_name = filename + ".pgm";
        if (m_n_channels == 3) {
                magic_head = "P6";
                extended_name = filename + ".ppm";
        }

        fout.open(extended_name.c_str(), ios::out | ios::binary);
        if (!fout.good()) {
                cerr << "Error opening file " << extended_name << " for output!" << endl;
                return false;
        }

        fout << magic_head << "\n";
        fout << m_width << " " << m_height << " 255\n";
        for (int y = 0; y < m_height; ++y) {
                const uchar *row_data = data(y);
                fout.write(reinterpret_cast<const char*>(row_data), m_width*m_n_channels*sizeof(uchar));
        }
        fout.close();

        return true;
}

bool Image::read_pnm(const std::string& filename)
{
        ifstream fin(filename.c_str(), ios::in | ios::binary);
        if (!fin.good()) {
                cerr << "Error opening PNM file " << filename << endl;
                return false;
        }

        int width;
        int height;
        int max_val;
        int n_channels = 1;
        string head = "00";
        head[0] = fin.get();
        head[1] = fin.get();
        if (head == "P5") {
                clog << "Loading PGM Binary" << endl;
                n_channels = 1;
        } else if (head == "P6") {
                clog << "Loading PPM Binary" << endl;
                n_channels = 3;
        } else {
                cerr << "File " << filename << " is not a Binary PGM or PPM!" << endl;
                return false;
        }

        fin >> width;
        fin >> height;
        fin >> max_val;
        if (fin.peek() == '\n')
                fin.get();

        int step = width * n_channels;
        uchar *new_data = new uchar[step*height];
        for (int y = 0; y < height; ++y) {
                fin.read(reinterpret_cast<char*>(new_data + y*step), step*sizeof(uchar));
        }
        fin.close();

        delete [] m_data;
        m_data = new_data;
        m_width = width;
        m_height = height;
        m_step = step;
        m_n_channels = n_channels;

        return true;
}

short *Image::deriv_x() const
{
        if (m_n_channels == 3) {
                cerr << "Image derivatives only implemented for grayscale images!" << endl;
                return nullptr;
        }

        short *dx = new short[m_width * m_height];
        for (int y = 0; y < m_height; ++y) {
                const uchar *row = this->data(y);
                short *drow = dx + y * m_width;
                drow[0] = 0;
                for (int x = 1; x < m_width - 1; ++x) {
                        drow[x] = row[x + 1] - row[x - 1];
                }
                drow[m_width - 1] = 0;
        }

        return dx;
}

short *Image::deriv_y() const
{
        if (m_n_channels == 3) {
                cerr << "Image derivatives only implemented for grayscale images!" << endl;
                return nullptr;
        }

        short *dy = new short[m_width * m_height];

        memset(dy, 0, m_width * sizeof(*dy));
        for (int y = 1; y < m_height - 1; ++y) {
                const uchar *rowm = this->data(y - 1);
                const uchar *rowp = this->data(y + 1);
                short *drow = dy + y * m_width;
                for (int x = 0; x < m_width; ++x) {
                        drow[x] = rowp[x] - rowm[x];
                }
        }
        memset(dy + (m_height - 1) * m_width, 0, m_width * sizeof(*dy));

        return dy;
}

void Image::rotate(Image *rotated, double theta, double tx, double ty) const
{
        if (m_n_channels != 1) {
                cerr << "Rotate only works on grayscale images!" << endl;
                return;
        }
        rotated->to_grayscale();

        double ct = cos(theta);
        double st = sin(theta);
        double tx_inv = -ct * tx + st * ty;
        double ty_inv = -st * tx - ct * ty;

        int wp = rotated->w();
        int hp = rotated->h();

        for (int yp = 0; yp < hp; ++yp) {
                uchar *row_p = rotated->data(yp);
                for (int xp = 0; xp < wp; ++xp) {
                        double x = ct * xp - st * yp + tx_inv;
                        double y = st * xp + ct * yp + ty_inv;

                        int x0 = static_cast<int>(x);
                        int y0 = static_cast<int>(y);

                        int value = 0;
                        if (x0 < 0 || y0 < 0 || x0 >= m_width || y0 >= m_height) {
                                value = 0;
                        } else {
                                const uchar *row = this->data(y0);
                                value = row[x0];
                        }

                        row_p[xp] = value;
                }
        }
}

void Image::rotate_centered(Image *rotated, double theta) const
{
        double ct = cos(theta);
        double st = sin(theta);
        double hw = m_width / 2.0;
        double hh = m_height / 2.0;
        double hwp = rotated->w() / 2.0;
        double hhp = rotated->h() / 2.0;

        double tx_cap = -ct * hw - st * hh + hwp;
        double ty_cap =  st * hw - ct * hh + hhp;
        this->rotate(rotated, theta, tx_cap, ty_cap);
}

void Image::smooth_x(float sigma)
{
        if (m_n_channels != 1) {
                cerr << "Smooth-x only works on grayscale images!" << endl;
                return;
        }

        int k = 0;
        unique_ptr<float []> kernel(gaussian_kernel(sigma, &k));

        int l = k / 2;
        unique_ptr<float []>  buffer(new float[m_width + 2 * l]);

        for (int y = 0; y < m_height - 1; ++y) {
                copy_to_buffer(buffer.get(), this->data(y), m_width, l, 1);
                convolve_buffer(buffer.get(), m_width, kernel.get(), k);
                copy_from_buffer(this->data(y), buffer.get(), m_width, 1);
        }
}

void Image::smooth_y(float sigma)
{
        if (m_n_channels != 1) {
                cerr << "Smooth-x only works on grayscale images!" << endl;
                return;
        }

        int k = 0;
        unique_ptr<float []> kernel(gaussian_kernel(sigma, &k));

        int l = k / 2;
        unique_ptr<float []>  buffer(new float[m_height + 2 * l]);

        for (int x = 0; x < m_width - 1; ++x) {
                copy_to_buffer(buffer.get(), m_data + x, m_height, l, m_step);
                convolve_buffer(buffer.get(), m_height, kernel.get(), k);
                copy_from_buffer(m_data + x, buffer.get(), m_height, m_step);
        }
}

void Image::smooth(float sigma_x, float sigma_y)
{
        smooth_x(sigma_x);
        smooth_y(sigma_y);
}

vector<Keypoint> Image::harris_corners(float threshold, float k, float sigma)
{       
        //(a) Computing x and y derivatives
        short *Ix = this->deriv_x();
        short *Iy = this->deriv_y(); 

        //(b) Computing squared values and x*y value of derivatives
        int length = m_width*m_height;
        short *Ix2 = new short[m_width * m_height];
        short *Iy2 = new short[m_width * m_height];
        short *Ixy = new short[m_width * m_height];
        for (int i = 0; i < length; ++i)
        {
              Ix2[i] = Ix[i]*Ix[i];//A values for M matrix
              Iy2[i] = Iy[i]*Iy[i];//C values for M matrix
              Ixy[i] = Ix[i]*Iy[i];//B values for M matrix               
        }

        //(c) Getting the image representations of calculated values and convolving them to get mini M matrices
        //consist of A B C
        Image *dx_img = short_to_image(Ix2, m_width, m_height);//A
        Image *dy_img = short_to_image(Iy2, m_width, m_height);//C
        Image *dxy_img = short_to_image(Ixy, m_width, m_height);//B
        //Smoothing them to get weighted images for derivatives Ix2 Iy2 Ixy
        dx_img->smooth(sigma,sigma);
        dy_img->smooth(sigma,sigma);
        dxy_img->smooth(sigma,sigma);

        //taking the pixel data arrays of derived and gaussian smoothed images 
        uchar* A = dx_img->data();//here we could put these values to a big array called M matrice as explained in harris corners formula but
        uchar* C = dy_img->data();//this would allocate more memory and for efficiency in the Score calculation we directly use these 3 arrays
        uchar* B = dxy_img->data();
        
        //(d) Creating an empty keypoint vector.
        vector<Keypoint> key_points;

        //e Calculating the score value for each pixel and store them in a float array (S)
        float* S = new float[m_width * m_height]();
        for (int y = 0; y < m_height; ++y)
        {
            for(int x = 0; x < m_width; ++x){
                int loc = y*m_width+x;
                float det = (A[loc]*C[loc]) - (B[loc]*B[loc]);
                float trace = A[loc] + C[loc];

                S[loc] = det - k * (trace * trace);//calculating score value for each pixel
            }
        }

       
        //f Check the Score values of each pixel according to Harris Cornerness controls.
        for (int y = 1; y < m_height-1; ++y)//we start from 1 and go to m_height-1 to avoid borders.
        {
            for (int x = 1; x < m_width-1; ++x)//we start from 1 and go to m_width-1 to avoid borders.
            {
                float current = S[y*m_width+x];
                float right = S[y*m_width+x+1];
                float left = S[y*m_width+x-1];
                float up = S[(y-1)*m_width+x];
                float down = S[(y+1)*m_width+x];
                if (current > left && 
                    current > right && 
                    current > up && 
                    current > down &&
                    current > threshold)//also compare with threshold
                {              
                       //if the condition is satisfied created a keypoint and add it to "points" vector
                        key_points.push_back ({(float) x,(float) y, current});                        
                }

            }
        }
        
        //converting image to rgb and putting green dots to the keypoints to visualize them      
        this->to_rgb();
        for (int i = 0; i < key_points.size(); ++i)
        {
            uchar *row_old = this->m_data + this->m_step * (int) key_points[i].y;
            row_old[(3* (int) key_points[i].x )+ 1] = 255;
        }

        delete [] Ix;
        delete [] Iy;
        delete [] Ix2;
        delete [] Iy2;
        delete [] Ixy;
        delete [] S;
        delete dx_img;
        delete dy_img;
        delete dxy_img;
        
        return key_points;  
}

}
