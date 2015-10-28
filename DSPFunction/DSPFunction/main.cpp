//
//  main.cpp
//  DSPFunction
//
//  Created by Linh Nguyen on 10/23/15.
//  Copyright (c) 2015 Linh Nguyen. All rights reserved.
//

#include <iostream>
#include "DSPFunction.h"
#include <random>


#include <string>
#include <math.h>


int main()
{
    //Testing functions for 1 second samples sample at 1Khz
    // Output should be total sum = 14 and PIR sum = 4 for XBand 7out of 10, PIR 4 out of 10
    //WalkingAwaySensors.txt first samples
    
    int lenPIR = 100;
    int lenXBand = 3000;
    
    float PIR[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0};
    float XBand[] = {1995,1989,1996,1998,1998,1996,1998,1998,1996,1998,1999,1998,1998,1999,1998,1997,1999,1998,1997,1999,1997,1996,1998,1998,1996,1998,1999,1997,1997,1999,1998,1997,1998,1998,1997,1998,1999,1996,1998,1998,1996,1996,1998,1997,1996,1997,1996,1996,1998,1998,1997,1997,1999,1999,1998,1999,2000,1998,1998,2000,1998,1998,1999,1998,1998,1998,1998,1996,1997,1997,1996,1997,1997,1995,1996,1998,1996,1996,1998,1997,1996,1998,1998,1997,1998,1998,1997,1998,1998,1997,1996,1997,1997,1995,1996,1997,1996,1996,1997,1997,1996,1999,1997,1998,1999,1997,1998,1999,1997,1999,1999,1998,1998,1998,1998,1999,1997,1998,1996,1997,1997,1998,1996,1997,1996,1996,1997,1996,1997,1997,1997,1997,1997,1999,1998,1996,1999,1998,1997,1998,1998,1996,1997,1998,1996,1995,1998,1996,1998,1998,1998,1997,1997,1997,2004,2004,2001,2002,2002,2004,2002,2004,2004,2003,2004,2004,2003,2004,2004,2004,2005,2005,2004,2006,2004,2005,2005,2004,2005,2005,2003,2004,2003,2002,2004,2003,2002,2004,2002,2003,2004,2003,2005,2004,2004,2005,2003,2005,2006,2004,2006,2006,2004,2005,2004,2003,2005,2004,2002,2003,2003,2002,2004,2004,2002,2004,2004,2003,2005,2005,2004,2006,2004,2005,2005,2005,2006,2004,2003,2005,2004,2003,2005,2006,2007,2006,2006,2007,2007,2006,2007,2008,2006,2007,2009,2008,2008,2008,2008,2009,2008,2008,2008,2008,2007,2008,2007,2007,2007,2008,2006,2007,2007,2006,2007,2008,2007,2007,2007,2006,2009,2007,2007,2008,2007,2008,2009,2007,2009,2008,2006,2008,2008,2006,2007,2007,2006,2006,2007,2005,2005,2006,2005,2006,2007,2006,2007,2007,2006,2007,2008,2006,2007,2007,2006,2007,2008,2006,2007,2006,2006,2006,2007,2005,2007,2006,2005,2007,2006,2007,2008,2006,2007,2008,2006,2007,2007,2006,2008,2007,2006,2008,2005,2006,2007,2005,2005,2006,2005,2005,2006,2006,2006,2006,2005,2007,2006,2007,2008,2006,2007,2007,2006,2008,2006,2006,2006,2004,2004,2005,2004,2005,2004,2005,2004,2005,2006,2004,2006,2006,2005,2006,2006,2004,2006,2007,2005,2006,2007,2005,2006,2007,2004,2004,2006,2004,2004,2005,2004,2004,2005,2004,2004,2005,2005,2005,2006,2004,2005,2006,2006,2005,2007,2005,2006,2006,2004,2006,2005,2004,2004,2005,2002,2004,2005,2004,2005,2005,2004,2005,2006,2005,2006,2006,2005,2006,2006,2005,2006,2006,2005,2005,2007,2004,2006,2005,2003,2005,2005,2003,2005,2005,2004,2005,2005,2004,2005,2005,2005,2006,2005,2004,2005,2005,2004,2005,2005,2005,2006,2004,2003,2005,2004,2004,1998,1997,1996,1998,1999,1997,1998,1999,1997,1999,1999,1998,2001,1998,1999,2000,1998,1999,1998,1997,1999,1998,1997,1998,1998,1997,1998,1998,1997,1998,1998,1997,1998,1999,1999,1997,1999,1999,1997,1999,1997,1998,1999,1998,1996,1998,1996,1997,1998,1997,1996,1997,1995,1995,1997,1996,1997,1997,1996,1996,1997,1997,1997,1998,1996,1997,1997,1997,1997,1997,1996,1996,1997,1995,1996,1997,1996,1995,1997,1995,1997,1997,1996,1997,1998,1996,1997,1998,1999,1996,1998,1995,1996,1998,1997,1997,1998,1995,1995,1997,1995,1995,1996,1996,1996,1996,1996,1994,1997,1996,1997,1998,1995,1997,1998,1997,1996,1998,1996,1996,1997,1995,1994,1997,1995,1995,1996,1995,1996,1997,1996,1995,1997,1995,1995,1997,1996,1997,1997,1996,1997,1996,1995,1997,1996,1995,1995,1995,1994,1995,1994,1993,1995,1994,1994,1994,1994,1994,1995,1995,1993,1995,1994,1995,1996,1994,1994,1994,1993,1993,1993,1992,1993,1992,1991,1992,1992,1991,1993,1993,1992,1992,1993,1992,1993,1993,1992,1993,1993,1992,1992,1993,1991,1993,1991,1991,1992,1991,1990,1991,1991,1990,1991,1990,1990,1991,1990,1990,1992,1990,1991,1992,1991,1991,1993,1991,1991,1992,1990,1989,1991,1982,1989,1984,1990,1991,1989,1990,1991,1990,1990,1992,1998,2000,2001,2000,2001,2001,2000,2002,2002,2001,2002,2001,2001,2002,2000,2001,2001,2000,2002,2001,2000,2002,2001,2001,2002,2001,2003,2002,2002,2003,2001,2003,2003,2002,2003,2004,2002,2004,2004,2002,2004,2003,2003,2004,2002,2002,2003,2001,2003,2003,2002,2003,2002,2003,2004,2002,2003,2004,2002,2004,2003,2003,2005,2003,2004,2004,2003,2004,2005,2003,2005,2004,2003,2005,2004,2005,2006,2005,2006,2006,2005,2007,2005,2005,2007,2005,2007,2007,2005,2005,2006,2005,2006,2006,2005,2006,2006,2005,2006,2007,2005,2006,2007,2006,2006,2008,2006,2006,2008,2007,2007,2010,2006,2007,2009,2007,2007,2009,2008,2007,2009,2008,2008,2009,2009,2008,2010,2009,2009,2011,2010,2009,2012,2011,2010,2012,2011,2010,2011,2012,2010,2011,2011,2011,2012,2011,2011,2012,2010,2012,2011,2012,2013,2011,2013,2013,2012,2014,2013,2013,2014,2013,2014,2015,2014,2014,2014,2013,2013,2014,2014,2015,2015,2013,2015,2015,2014,2018,2015,2014,2017,2016,2016,2018,2016,2015,2017,2016,2015,2017,2016,2017,2015,2016,2018,2015,2017,2016,2015,2017,2017,2017,2018,2018,2018,2018,2018,2018,2018,2018,2019,2018,2018,2020,2018,2018,2019,2018,2018,2020,2018,2019,2021,2018,2020,2020,2019,2021,2019,2020,2021,2019,2000,2003,2004,2003,2003,2004,2002,2003,2003,2002,2003,2001,2002,2003,2001,2003,2002,2001,2003,2003,2001,2004,2003,2003,2005,2004,2003,2005,2003,2004,2005,2003,2005,2004,2003,2004,2002,2002,2004,2002,2002,2003,2003,2003,2003,2003,2004,2002,2004,2003,2003,2005,2004,2004,2003,2003,2006,2003,2003,2005,2002,2003,2004,2003,2004,2002,2004,2004,2003,2005,2005,2003,2006,2003,2004,2006,2005,2006,2006,2004,2006,2004,2004,2006,2003,2004,2004,2004,2005,2003,2003,2003,2002,2004,2004,2003,2005,2004,2003,2005,2004,2004,2005,2004,2007,2004,2004,2005,2003,2005,2003,2003,2003,2002,2004,2005,2002,2005,2002,2004,2006,2005,2005,2006,2003,2005,2006,2004,2006,2005,2007,2005,2004,2006,2004,2003,2005,2003,2004,2005,2002,2003,2004,2002,2004,2004,2003,2005,2005,2005,2005,2003,2005,2005,2005,2005,2006,2005,2005,2003,2005,2004,2002,2004,2002,2003,2004,2003,2003,2005,2003,2004,2005,2003,2005,2005,2003,2005,2005,2003,2004,2005,2003,2006,2003,2003,2005,2002,2003,2003,2002,2004,2004,2002,2003,2004,2002,2005,2004,2003,2005,2004,2005,2005,2005,2004,2005,2004,2004,2004,2002,2003,2003,2002,2004,2003,2002,2004,2002,2003,2004,2003,2003,2004,2004,2003,2005,2004,2003,2005,2003,2004,1978,1980,1979,1977,1978,1978,1978,1979,1977,1980,1978,1978,1979,1978,1978,1979,1979,1978,1980,1980,1979,1980,1979,1980,1981,1979,1979,1979,1978,1978,1979,1978,1976,1978,1976,1976,1977,1976,1977,1977,1976,1978,1977,1977,1979,1978,1978,1979,1977,1978,1979,1976,1977,1977,1975,1977,1977,1976,1977,1976,1975,1976,1977,1976,1977,1977,1976,1977,1977,1977,1978,1977,1977,1977,1976,1978,1977,1976,1976,1976,1975,1976,1974,1974,1975,1974,1974,1974,1974,1976,1974,1976,1975,1975,1975,1974,1976,1974,1976,1975,1975,1975,1974,1975,1974,1974,1974,1972,1974,1973,1973,1974,1972,1974,1972,1974,1974,1973,1974,1974,1976,1975,1974,1975,1974,1974,1975,1973,1974,1973,1972,1973,1972,1971,1972,1971,1970,1972,1971,1972,1973,1972,1973,1974,1971,1974,1973,1972,1974,1973,1972,1974,1972,1973,1971,1970,1971,1970,1970,1972,1970,1970,1972,1970,1971,1972,1970,1972,1972,1971,1973,1973,1974,1974,1972,1974,1974,1971,1973,1972,1971,1972,1971,1970,1971,1969,1971,1970,1970,1971,1971,1972,1973,1971,1972,1973,1971,1972,1973,1971,1972,1972,1971,1971,1969,1971,1972,1969,1971,1971,1970,1970,1970,1969,1970,1971,1969,1971,1972,1971,1972,1973,1972,1974,1973,1971,1973,1972,1971,1972,1971,1970,2014,2015,2014,2013,2015,2015,2014,2015,2017,2015,2015,2018,2016,2016,2018,2017,2017,2019,2018,2017,2018,2018,2017,2017,2018,2018,2015,2017,2016,2015,2017,2017,2016,2017,2017,2016,2013,2016,2019,2017,2018,2019,2018,2020,2018,2019,2020,2018,2019,2020,2017,2018,2017,2017,2018,2017,2017,2019,2018,2018,2019,2019,2018,2021,2018,2020,2019,2019,2021,2020,2021,2021,2020,2020,2021,2020,2019,2018,2020,2020,2019,2018,2019,2018,2020,2019,2018,2020,2018,2020,2020,2020,2021,2021,2021,2022,2020,2021,2020,2019,2020,2020,2019,2020,2019,2020,2020,2019,2020,2022,2019,2021,2020,2020,2021,2023,2021,2023,2022,2021,2023,2022,2020,2022,2020,2020,2022,2020,2020,2019,2019,2019,2021,2020,2021,2021,2021,2020,2021,2022,2021,2022,2023,2021,2022,2021,2021,2022,2021,2021,2021,2020,2019,2021,2020,2021,2022,2021,2021,2022,2021,2022,2022,2022,2021,2023,2023,2022,2023,2023,2022,2023,2022,2022,2023,2021,2021,2021,2021,2021,2021,2020,2021,2022,2021,2021,2022,2022,2022,2023,2022,2023,2024,2023,2024,2023,2023,2024,2022,2023,2021,2022,2023,2022,2020,2023,2021,2022,2024,2022,2023,2025,2022,2023,2026,2024,2025,2027,2024,2025,2025,2024,2025,2026,2024,2024,2025,2023,2023,2025,2023,2023,2033,2031,2033,2033,2033,2034,2033,2035,2034,2033,2034,2035,2033,2034,2035,2033,2035,2035,2032,2034,2032,2032,2033,2032,2031,2032,2031,2032,2033,2032,2033,2033,2031,2032,2033,2032,2033,2033,2031,2032,2033,2031,2032,2032,2032,2030,2032,2031,2030,2032,2030,2030,2032,2031,2030,2031,2030,2030,2031,2031,2030,2031,2031,2030,2031,2030,2029,2031,2031,2028,2030,2030,2028,2028,2029,2027,2028,2029,2027,2027,2029,2028,2028,2029,2029,2028,2029,2029,2028,2028,2029,2028,2027,2028,2028,2026,2027,2026,2025,2027,2026,2025,2028,2025,2026,2027,2025,2027,2026,2026,2028,2027,2025,2027,2028,2026,2024,2026,2024,2024,2026,2025,2024,2023,2023,2023,2023,2024,2024,2023,2024,2024,2024,2025,2024,2022,2025,2025,2022,2022,2025,2022,2023,2022,2022,2021,2021,2021,2021,2020,2022,2020,2019,2021,2020,2020,2021,2021,2020,2022,2022,2020,2022,2021,2020,2021,2021,2019,2021,2019,2019,2020,2019,2018,2019,2019,2018,2018,2018,2018,2017,2019,2019,2018,2020,2019,2018,2019,2020,2018,2018,2019,2017,2016,2018,2015,2017,2017,2016,2017,2017,2015,2016,2017,2015,2017,2017,2015,2016,2017,2016,2016,2017,2015,2016,2016,2015,2015,2015,2014,2015,2015,2014,2013,2015,2013,2013,2014,2013,2013,2013,2012,1987,1987,1987,1987,1986,1989,1988,1986,1988,1986,1987,1987,1987,1988,1987,1987,1988,1986,1988,1987,1987,1987,1988,1986,1988,1987,1988,1989,1987,1989,1988,1987,1989,1988,1987,1989,1987,1988,1989,1988,1987,1989,1988,1988,1989,1987,1989,1989,1988,1990,1989,1988,1989,1990,1989,1990,1991,1989,1990,1991,1989,1989,1991,1989,1989,1990,1989,1989,1990,1989,1990,1991,1988,1989,1990,1990,1989,1991,1989,1990,1991,1990,1991,1993,1990,1991,1991,1990,1991,1991,1990,1991,1991,1991,1992,1991,1991,1993,1992,1992,1992,1991,1992,1993,1992,1993,1993,1991,1992,1994,1992,1993,1992,1992,1992,1994,1991,1992,1993,1992,1993,1993,1991,1992,1993,1991,1993,1993,1991,1994,1992,1992,1993,1993,1993,1994,1992,1994,1993,1992,1994,1993,1993,1994,1993,1992,1993,1993,1994,1994,1994,1992,1995,1994,1994,1996,1994,1993,1995,1993,1995,1996,1994,1995,1996,1994,1996,1995,1995,1994,1994,1994,1997,1994,1995,1996,1995,1996,1996,1994,1997,1997,1997,1996,1995,1997,1996,1996,1995,1996,1995,1994,1996,1995,1994,1997,1996,1998,1997,1997,1997,1998,1996,1997,1998,1997,1998,1998,1997,1997,1999,1997,1997,1997,1997,1998,1998,1996,1997,1998,1996,1998,1997,1997,1997,1997,1997,1998,1998,1999,2000,1998,2030,2019,2032,2031,2030,2030,2030,2030,2031,2029,2028,2032,2030,2031,2030,2032,2029,2030,2031,2029,2030,2031,2029,2031,2032,2030,2032,2030,2032,2032,2030,2032,2032,2030,2031,2032,2030,2031,2030,2030,2031,2029,2030,2030,2030,2031,2030,2030,2032,2029,2030,2031,2030,2031,2031,2031,2031,2030,2029,2031,2028,2030,2030,2028,2029,2029,2027,2028,2028,2028,2029,2028,2030,2029,2028,2030,2028,2028,2030,2029,2029,2029,2028,2029,2027,2028,2028,2027,2028,2028,2027,2028,2027,2027,2027,2026,2028,2027,2027,2028,2027,2027,2028,2026,2027,2027,2027,2026,2024,2026,2026,2025,2026,2025,2024,2025,2023,2024,2025,2023,2024,2025,2024,2025,2025,2023,2025,2024,2024,2025,2023,2023,2025,2023,2023,2023,2022,2022,2022,2021,2022,2021,2022,2023,2021,2022,2023,2021,2023,2023,2020,2024,2022,2021,2023,2022,2020,2019,2021,2020,2019,2021,2018,2021,2019,2019,2019,2018,2020,2018,2020,2018,2021,2019,2019,2020,2018,2020,2018,2019,2020,2017,2018,2019,2015,2017,2016,2015,2016,2016,2014,2016,2014,2015,2016,2015,2015,2017,2014,2016,2016,2014,2015,2016,2014,2015,2015,2014,2015,2014,2012,2013,2012,2011,2012,2012,2011,2012,2010,2012,2012,2010,2011,2011,2010,2012,2010,2010,2011,2010,2010,2010,1928,1930,1929,1926,1927,1925,1924,1925,1923,1923,1923,1922,1921,1921,1920,1920,1921,1919,1920,1921,1918,1919,1918,1918,1917,1916,1916,1914,1912,1913,1911,1911,1911,1909,1910,1908,1910,1909,1906,1909,1908,1906,1907,1906,1906,1906,1905,1905,1905,1903,1903,1902,1901,1902,1899,1900,1898,1898,1897,1895,1898,1896,1894,1895,1895,1894,1895,1893,1893,1894,1892,1893,1893,1891,1893,1889,1889,1890,1887,1887,1888,1886,1886,1885,1884,1885,1885,1884,1885,1884,1883,1884,1884,1883,1884,1882,1883,1883,1881,1882,1882,1878,1880,1879,1876,1877,1877,1875,1876,1874,1874,1874,1874,1875,1873,1874,1874,1873,1874,1872,1874,1872,1872,1872,1871,1871,1870,1870,1869,1868,1868,1867,1868,1868,1866,1868,1866,1867,1866,1865,1866,1867,1864,1866,1866,1866,1867,1865,1864,1865,1864,1865,1865,1863,1863,1861,1862,1862,1860,1861,1858,1858,1860,1858,1860,1859,1858,1860,1859,1859,1859,1860,1859,1859,1859,1858,1857,1858,1858,1855,1855,1855,1854,1855,1855,1854,1855,1855,1855,1856,1856,1855,1857,1854,1856,1856,1855,1856,1856,1855,1856,1855,1855,1855,1853,1854,1854,1853,1853,1853,1852,1851,1853,1853,1853,1853,1854,1855,1853,1855,1853,1854,1855,1854,1854,1855,1853,1853,1853,1853,1853,1854,1893,1894,1893,1893,1893,1893,1894,1894,1894,1893,1897,1896,1896,1895,1898,1897,1896,1899,1898,1898,1896,1899,1897,1896,1898,1895,1896,1898,1894,1897,1899,1897,1900,1900,1899,1901,1901,1901,1903,1902,1902,1904,1903,1902,1905,1902,1902,1904,1902,1902,1903,1901,1900,1903,1901,1902,1902,1903,1903,1903,1904,1905,1905,1906,1906,1905,1908,1906,1906,1909,1907,1906,1908,1906,1905,1906,1907,1905,1907,1907,1906,1909,1908,1908,1911,1909,1910,1912,1912,1912,1913,1913,1912,1915,1914,1913,1915,1913,1912,1913,1912,1911,1912,1911,1911,1913,1911,1913,1914,1914,1914,1917,1915,1917,1918,1917,1919,1918,1918,1919,1918,1918,1919,1917,1917,1918,1916,1917,1918,1918,1919,1920,1919,1920,1922,1921,1922,1924,1923,1925,1926,1923,1926,1926,1925,1927,1926,1925,1926,1926,1924,1925,1924,1924,1925,1925,1926,1926,1926,1927,1928,1927,1929,1929,1928,1931,1931,1930,1931,1931,1931,1930,1931,1929,1930,1930,1929,1931,1929,1932,1932,1931,1933,1933,1934,1935,1935,1935,1938,1936,1938,1938,1938,1940,1940,1937,1940,1939,1938,1939,1937,1938,1938,1937,1938,1938,1937,1940,1938,1939,1941,1941,1942,1943,1942,1943,1945,1943,1944,1944,1944,1945,1942,1943,1943,1942,1942,1944,1943,1944,1945,1944,2016,2016,2017,2016,2018,2019,2018,2020,2019,2019,2020,2021,2019,2020,2021,2018,2019,2020,2016,2016,2017,2016,2017,2017,2018,2018,2018,2018,2020,2019,2019,2020,2021,2021,2021,2021,2021,2022,2021,2021,2022,2021,2020,2020,2017,2019,2019,2017,2018,2017,2017,2018,2018,2017,2019,2019,2018,2021,2020,2019,2020,2021,2020,2021,2021,2019,2020,2020,2017,2019,2018,2016,2018,2019,2016,2018,2020,2019,2018,2019,2021,2021,2020,2022,2021,2021,2023,2021,2021,2022,2020,2022,2020,2021,2020,2017,2019,2020,2017,2018,2018,2018,2018,2019,2018,2021,2020,2018,2018,2019,2020,2021,2020,2019,2021,2020,2019,2020,2019,2017,2018,2018,2017,2018,2018,2018,2018,2020,2019,2019,2020,2018,2019,2022,2021,2022,2022,2022,2022,2023,2020,2021,2022,2019,2020,2020,2017,2018,2017,2017,2019,2018,2017,2019,2019,2018,2020,2020,2020,2021,2019,2021,2021,2020,2021,2019,2019,2020,2018,2017,2017,2017,2016,2018,2018,2015,2019,2017,2018,2021,2018,2019,2023,2020,2021,2022,2021,2021,2021,2021,2022,2021,2019,2019,2019,2019,2018,2017,2017,2018,2018,2017,2018,2017,2018,2019,2018,2020,2020,2019,2021,2021,2019,2020,2020,2018};
//    float XBand[] = {2067,2065,2059,2056,2049,2046,2041,2038,2034,2029,2026,2023,2019,2016,2012,2011,1998,1993,1994,1996,1996,1994,1992,1990,1989,1990,1987,1985,1983,1981,1981,1980,1979,1976,1977,1976,1976,1976,1974,1973,1973,1974,1973,1974,1975,1974,1972,1973,1974,1975,1975,1976,1976,1975,1978,1979,1979,1978,1980,1982,1983,1983,1984,1985,1987,1989,1991,1990,1992,1993,1996,1996,1999,2001,2002,2004,2005,2007,2009,2012,2012,2013,2016,2018,2019,2021,2025,2027,2028,2031,2034,2036,2037,2039,2042,2044,2047,2049,2050,2053,2055,2058,2060,2061,2064,2066,2069,2071,2072,2073,2075,2079,2082,2083,2084,2085,2086,2088,2089,2092,2094,2095,2095,2099,2102,2102,2103,2105,2105,2108,2108,2108,2109,2111,2111,2103,2105,2107,2112,2112,2114,2115,2114,2115,2117,2116,2116,2114,2114,2115,2115,2113,2113,2113,2113,2112,2109,2109,2109,2106,2105,2104,2102,2099,2098,2098,2095,2091,2089,2089,2086,2081,2081,2078,2074,2070,2068,2066,2062,2058,2055,2053,2048,2045,2042,2037,2031,2030,2025,2019,2017,2014,2009,2003,1999,1996,1992,1985,1983,1978,1974,1967,1960,1961,1953,1949,1945,1939,1935,1929,1925,1920,1913,1911,1901,1896,1893,1890,1885,1881,1878,1873,1867,1862,1859,1854,1850,1845,1839,1836,1833,1828,1823,1818,1816,1813,1808,1805,1800,1797,1794,1789,1787,1783,1779,1777,1775,1772,1768,1767,1762,1761,1759,1756,1756,1752,1751,1750,1748,1747,1747,1744,1745,1744,1741,1742,1742,1741,1740,1742,1741,1742,1743,1740,1733,1739,1741,1744,1747,1749,1750,1752,1755,1756,1758,1760,1762,1764,1767,1769,1770,1773,1776,1779,1780,1784,1787,1789,1792,1796,1797,1801,1805,1806,1809,1812,1816,1818,1822,1826,1827,1831,1835,1836,1840,1844,1847,1847,1851,1855,1857,1859,1862,1865,1869,1870,1872,1876,1880,1881,1882,1884,1888,1891,1892,1895,1898,1901,1901,1903,1907,1909,1909,1912,1915,1916,1917,1920,1923,1924,1924,1927,1930,1932,1933,1933,1935,1937,1939,1941,1942,1942,1943,1946,1948,1948,1949,1951,1953,1955,1957,1957,1957,1959,1959,1962,1961,1962,1963,1965,1967,1967,1967,1968,1968,1970,1971,1971,1971,1971,1974,1975,1974,1973,1975,1976,1977,1977,1976,1978,1979,1978,1977,1980,1980,1981,1981,1980,1980,1981,1982,1982,1982,1982,1981,1982,1983,1983,1983,1982,1983,1984,1984,1985,1984,1983,1984,1985,1985,1984,1984,1984,1985,1985,1985,1987,1986,1984,1985,1986,1987,1986,1985,1986,1987,1986,1986,1987,1988,1986,1986,1988,1989,1986,1986,1988,1988,1987,1989,1989,1987,1988,1989,1988,1989,1990,1988,1990,1990,1988,1990,1990,1990,1990,1991,1990,1990,1992,1992,1991,1991,1992,1993,1993,1993,1992,1992,1993,1994,1993,1992,1993,1994,1994,1992,1995,1993,1994,1995,1995,1994,1996,1993,1996,1996,1994,1995,1996,1994,1994,1994,1996,1994,1994,1996,1995,1994,1995,1996,1995,1994,1996,1994,1994,1996,1996,1994,1996,1997,1997,2003,2009,2013,2014,2016,2019,2018,2017,2016,2018,2018,2016,2016,2019,2016,2016,2017,2017,2015,2017,2017,2015,2016,2017,2015,2017,2017,2015,2016,2017,2015,2016,2017,2015,2016,2018,2017,2016,2017,2016,2017,2015,2017,2017,2018,2018,2017,2019,2017,2020,2018,2020,2019,2020,2021,2019,2020,2021,2020,2019,2021,2022,2021,2022,2022,2022,2022,2023,2022,2023,2024,2024,2024,2025,2025,2026,2026,2026,2026,2028,2027,2018,2017,2021,2023,2026,2028,2028,2029,2030,2030,2030,2030,2031,2032,2030,2031,2032,2032,2031,2031,2031,2031,2032,2031,2031,2033,2032,2031,2031,2032,2031,2032,2032,2030,2031,2032,2030,2030,2030,2029,2017,2021,2024,2026,2025,2027,2028,2027,2026,2028,2028,2026,2026,2026,2026,2024,2025,2026,2023,2024,2025,2023,2022,2024,2022,2020,2020,2021,2019,2018,2019,2018,2016,2016,2016,2014,2015,2013,2012,1982,1950,1972,1988,1998,2001,2004,2007,1980,1956,1977,1989,1995,2000,2003,2003,2004,2005,2003,2002,2002,2000,2000,1999,1999,1997,1995,1995,1994,1993,1992,1991,1990,1987,1988,1985,1983,1983,1984,1981,1980,1979,1976,1976,1977,1975,1973,1973,1971,1969,1970,1966,1967,1966,1963,1962,1962,1961,1959,1958,1958,1956,1954,1954,1953,1950,1949,1950,1947,1945,1945,1944,1942,1941,1941,1940,1938,1937,1937,1937,1934,1934,1934,1933,1930,1929,1930,1929,1927,1928,1925,1927,1925,1924,1925,1924,1922,1921,1921,1921,1920,1919,1920,1918,1918,1916,1917,1916,1915,1915,1915,1916,1913,1913,1913,1914,1912,1911,1911,1911,1911,1911,1910,1910,1910,1910,1910,1908,1908,1909,1908,1906,1908,1907,1905,1906,1908,1906,1905,1905,1906,1905,1905,1904,1906,1906,1904,1903,1866,1872,1882,1890,1896,1901,1902,1902,1903,1903,1905,1904,1902,1903,1902,1904,1903,1902,1901,1901,1901,1901,1899,1898,1900,1900,1898,1898,1898,1898,1896,1896,1896,1895,1895,1897,1896,1894,1894,1895,1893,1893,1893,1894,1891,1891,1892,1892,1890,1889,1889,1890,1888,1888,1889,1888,1887,1887,1887,1888,1886,1886,1885,1886,1886,1885,1884,1886,1886,1884,1883,1884,1884,1884,1884,1882,1883,1884,1884,1882,1882,1882,1881,1882,1880,1880,1881,1881,1880,1879,1878,1879,1879,1879,1879,1877,1876,1870,1865,1861,1858,1857,1858,1857,1856,1855,1856,1857,1858,1857,1857,1857,1857,1857,1857,1857,1857,1859,1856,1857,1858,1859,1858,1858,1859,1860,1858,1859,1860,1860,1858,1859,1860,1860,1858,1858,1860,1860,1860,1859,1860,1861,1859,1861,1859,1860,1861,1860,1859,1860,1860,1860,1859,1859,1860,1861,1859,1859,1859,1861,1860,1859};
//    
    Detection detect;
    float noise_Var = 5;
    float signal_Var = 6.8273e5;
    float Vt = 4.4412;
    int Fs = 3000;
    int detector;
    int PIRsum;
    detect.SetValues(XBand, Fs);
    int mean = detect.Mean();
    detector = detect.Detector(XBand, PIR, lenXBand, lenPIR, noise_Var, Vt, signal_Var,Fs,mean,&PIRsum);
    printf("%d\t",detector);
    printf("%d\t",PIRsum);
//    int Fs = 1000;
//    int Fo = 100;
//    int Fo1 = 10;
//    int *BinaryIntegrationXB;
//    int *BinaryIntegrationPIR;
//    int lendetection = 20;
//    int *detection;
//    int *stage1Detection;
//    BinaryIntegrationXB = (int*) calloc(Fo1, sizeof(int));
//    BinaryIntegrationPIR = (int*) calloc(Fo1, sizeof(int));
//    detection = (int*) calloc(lendetection, sizeof(int));
//    stage1Detection = (int*) calloc(Fo1, sizeof(int));
//    float *XFiltered;
//
//    
//    float noise_Var = 5; //nose_Var and Vt change accoordingly
//    float signal_Var = 6.8273e5;
//    float Vt = 4.4412; //recalculate Vt when noise_Var change
//    float lrt =(2*pow(noise_Var,2)*pow(signal_Var,2)/(pow(signal_Var,2) - pow(noise_Var,2)))*log(Vt) - (Fs/Fo)*log(noise_Var/signal_Var);
//    
//    detect.SetValues(XBand, sizeof(XBand)/sizeof(XBand[0]));
//    XFiltered = detect.WMA(XBand,lenXBand,9);
//    
//
//    for (int observation = 0;observation < Fo; observation++)
//    {
//        float sumxobssq = 0;
//        for (int i = 0; i < Fs/Fo; i++)
//        {
//            sumxobssq = sumxobssq + pow(XFiltered[i + (Fs/Fo)*observation],2);
//        }
//        if (sumxobssq > lrt)
//        {
//            stage1Detection[observation] = 1;
//        }
//        else
//            stage1Detection[observation] = 0;
//    }
//
//
//    for (int observation = 0; observation < Fo1;observation++)
//    {
//        float sumXB = 0;
//        float sumPIR = 0;
//        for (int i = 0; i < Fo1; i++)
//        {
//            sumXB = sumXB + stage1Detection[i+observation*Fo1];
//            sumPIR = sumPIR + PIR[i+observation*Fo1];
//
//        }
//        if (sumXB > 7)
//            BinaryIntegrationXB[observation] = 1;
//        else
//            BinaryIntegrationXB[observation] = 0;
//        if (sumPIR > 7)
//            BinaryIntegrationPIR[observation] = 1;
//        else
//            BinaryIntegrationPIR[observation] = 0;
//    }
//    
//    for (int i = 0; i < lendetection/2; i++)
//    {
//        detection[i] = BinaryIntegrationXB[i];
//        detection[lendetection-i-1] = BinaryIntegrationPIR[i];
//    }
    
    
    // print out
//    for (int i = 0; i < 20;i++)
//    {
//        printf("%d \t",detector[i]);
//    }
//    printf("\n");

}


//    evariance = detec.CalculateSampleVariane();
//    double sampledevi = detec.GetSampleStandardDeviation();
//    double devi = detec.GetStandardDeviation();
//    float A[] = {1,2,3,4};
//    float B[] = {2,3,4,5};
//    float *C;
//    int lenC;
//    
//    Detection detec;
//    
//    C = detec.conv(A,B,4,4,&lenC);
//    
//    for (int i = 0; i < lenC;i++)
//    {
//        printf("%5.3f \t",C[i]);
//        printf("\n\n");
//    }
//    
////
//    double arrNumbers[] =
//    {
//        15.17, 16.94, 14.94, 14.99, 13.77, 13.75,
//        12.67, 12.14, 12.59, 12.48, 14.81, 14.29,
//        12.74, 12.52, 11.65, 12.24, 11.42, 12.25,
//        12.72, 11.64, 11.09, 11.22, 11.50, 11.36,
//        11.84, 12.18, 11.04, 10.90, 11.80, 11.84,
//    };
//    
//    
//    detec.SetValues(arrNumbers, sizeof(arrNumbers) / sizeof(arrNumbers[0]));
//    
//    double mean = detec.CalculateMean();
//    double variance = detec.CalculateVariane();
//    double sampl
//    char buf[1024];
//    sprintf(buf, "Total Numbers\t\t\t: %10lu\n", sizeof(arrNumbers) / sizeof(arrNumbers[0]));
//    std::cout << buf;
//    sprintf(buf, "Mean\t\t\t\t: %10.5lf\nPopulation Variance\t\t: %10.4lf\n", mean, variance);
//    std::cout << buf;
//    sprintf(buf, "Sample variance\t\t\t: %10.4lf\n", samplevariance);
//    std::cout << buf;
//    sprintf(buf, "Population Standard Deviation\t: %10.4lf\n", devi);
//    std::cout << buf;
//    sprintf(buf, "Sample Standard Deviation\t: %10.4lf\n", sampledevi);
//    std::cout << buf;
//
//    return 0;
//}


