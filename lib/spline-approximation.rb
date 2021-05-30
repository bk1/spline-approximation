#
# spline-approximation.rb
#
# (C) Karl Brodowsky (IT Sky Consulting GmbH) 2021
#
# Author:    bk1 (Karl Brodowsky)
#

require "rational"
require "long-decimal"
require "bigdecimal"
require "complex"

class SplineApproximation

  # find approximation g(x) for function on interval [x_min, x_max] with
  # x_max > x_min
  # n > 0 and integral number
  # h = (x_max - x_min)/n
  # x_i = x_min + i*h for i=0..n
  def initialize(x_min, x_max, n)
    if (x_min.kind_of? Complex) 
      raise "x_min must not be a complex number"
    end
    if (x_max.kind_of Complex)
      raise "x_max must not be a complex number"
    end
    if (x_min >= x_max)
      raise "x_min=#{x_min} must not be >= x_max=#{x_max}"
    end
    if (! (n.kind_of? Integer) || n <= 0)
      raise "n=#{n} must be an integer > 0"
    end
    @x_min = x_min
    @x_max = x_max
    # number of sub intervals
    @n = n
    # length of sub intervals
    @h = (x_max - x_min)/n
    @x = (0..n).map do |i|
      # calculate x_i
      x_min + i*h
    end
    @x_i = (-3..n+3).map do |i|
      x_min + i*h
    end      

    # rows 1..n with 
    @lin_equ = ((-1)..(n+1)).map do |i|
      Array.new(n+4, 0)
    end

    @count = 0
    
  end # initialize

  attr_reader :x_min, :x_max, :n, :h, :x, :count

  def x_i(i)
    @x_i[i+3]
  end

  def lin_equ(i, j)
    @lin_equ[i+1][j+1]
  end

  def f(x)
    if (x <= -2)
      0
    elsif (x <= -1)
      (x+2)**3
    elsif (x <= 0)
      xx = x+1
      ((-3*xx+3)*xx + 3)*xx +1
    elsif (x < 1)
      xx = 1-x
      ((-3*xx+3)*xx + 3)*xx +1
    elsif (x < 2)
      (2-x)**3
    else
      0
    end

  end

  def f_i(x,i)
    f((x-x_i(i))/h)
  end

  def aggregate(xi, eta, m = 1)
    (-1..n+1).each do |k|
      m_f_k_xi = f_i(xi, k)
      @lin_equ[k+1][n+3] += m_f_k_xi*eta
      (-1..n+1).each do |i|
        f_i_xi = f_i(xi, i)
        @lin_equ[k+1][i+1] += m_f_k_xi * f_i_xi
      end
    end
    @count += m
  end

  def calc_cub
    @cub_sums = Array.new(n+3, 0)
    i = -1
    @cub = @lin_equ.map do |line|
      k = -1
      cline = line.map do |element|
        c = element.abs**3
        if (k <= n+1)
          cub_sum[i+1]+=c
        end
        k += 1
        c
      end
      i += 1
      cline
    end
  end
  
  def solve
    calc_cub()
    (-1..n+1).each do |k|
      pivot_i = -1
      pivot_crit = 0
      (k..n+1).each do |i|
        pc = @cub[i+1,k+1] / @cub_sums[i+1]
        if (pc > pivot_crit)
          pivot_i = i
          pivot_crit = pc
        end
      end
      if (pivot_i != k)
        line = @lin_equ[i+1]
        @lin_equ[i+1] = @lin_equ[k+1]
        @lin_equ[k+1] = line
      end
      (k+1..n+1).each do |i|
        pivot = @lin_equ[k+1, k+1]
        if (@lin_equ[i+1, k+1] != 0)
          q = @lin_equ[i+1, k+1]/pivot
          @lin_equ[i+1, k+1] = 0
          (k+2..n+2).each do |l|
            @lin_equ[i+1, l+1] -= q*@lin_equ[k+1, l+1]
          end
        end
      end
      # inefficient but ok as first step
      calc_cub
    end
    
  end
  
end
