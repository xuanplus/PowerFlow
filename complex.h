#pragma once
#include <cmath>
#include <sstream>
#include <iostream>


// ������
class complex
{
public:

	float real;
	float imag;

	// ���캯��
	complex(float real, float imag);

	// ģֵ
	float modulus() const
	{
		return sqrtf(powf(this->real, 2) + powf(this->imag, 2));
	}

	// ����
	complex conj() const {
		return complex(real, -imag);
	}

	// ���ظ��������
	complex operator+(const complex& c);
	complex operator-(const complex& c);
	complex operator*(const complex& c);
	complex operator/(const complex& c);
	complex operator+=(const complex& c) { return *this = *this + c; }
	complex operator-=(const complex& c) { return *this = *this - c; }
	complex operator=(const complex& c) { this->real = c.real; this->imag = c.imag; return *this; }

	// ����==�����
	bool operator==(const complex& c) const {
		return this->real == c.real && this->imag == c.imag;
	}

	// ������������
	friend std::ostream& operator<<(std::ostream& os, const complex& c)
	{
		std::ostringstream temp;
		if (c.real == 0 && c.imag == 0)
		{
			temp << 0;
		}
		else
		{
			if (c.real != 0)
			{
				temp << c.real;
			}
			if (c.imag != 0)
			{
				temp << (c.imag > 0 ? "+j" : "-j") << abs(c.imag);
			}
		}
		os << temp.str();
		return os;
	}
};

complex::complex(float real, float imag)
{
	this->real = real;
	this->imag = imag;
}

inline complex complex::operator+(const complex& c)
{
	return complex(this->real + c.real, this->imag + c.imag);
}

inline complex complex::operator-(const complex& c)
{
	return complex(this->real - c.real, this->imag - c.imag);
}

inline complex complex::operator*(const complex& c)
{
	return complex(this->real * c.real - this->imag * c.imag, this->real * c.imag + this->imag * c.real);
}

inline complex complex::operator/(const complex& c)
{
	float denominator = c.real * c.real + c.imag * c.imag;
	return complex((this->real * c.real + this->imag * c.imag) / denominator,
		(this->imag * c.real - this->real * c.imag) / denominator);
}
