# coding: utf-8
#
# spline-approximation.rb
#
# (C) Karl Brodowsky (IT Sky Consulting GmbH) 2021
#
# Author:    bk1 (Karl Brodowsky)
#

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
    if (x_max.kind_of? Complex)
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
    len = x_max - x_min
    if (len.kind_of? Integer)
      len = len.to_f
    end
    @h = len/n
    @x = (0..n).map do |i|
      # calculate x_i
      x_min + i*h
    end
    @x_i = (-3..n+3).map do |i|
      x_min + i*h
    end
    @x_minus_3 = @x_i[0]
    @x_n_plus_3 = @x_i[n+6]

    # rows 1..n with 
    @lin_equ = ((-1)..(n+1)).map do |i|
      Array.new(n+4, 0)
    end

    @count = 0
    
  end # initialize

  attr_reader :x_min, :x_max, :n, :h, :x, :count, :results, :coefficients

  def to_s
    "SplineApproximation(x_min=#{x_min} x_max=#{x_max} n=#{n} count=#{count})"
  end

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
    # puts("f_i(x=#{x} i=#{i}) x_i=#{x_i(i)} h=#{@h}")
    f((x-x_i(i))/@h)
  end

  def aggregate_tolerant(xi, eta, m = 1)
    if (xi < x_min || xi > x_max)
      return
    end
    aggregate(xi, eta, m)
  end
    
  def aggregate(xi, eta, m = 1)
    if (xi < x_min || xi > x_max)
      raise ("x=#{x} is out of range [x_min, x_max]=[#{x_min}, #{x_max}]")
    end
    (-1..n+1).each do |k|
      m_f_k_xi = m*f_i(xi, k)
      # puts("#{m}*f_#{k}(#{xi})=#{m_f_k_xi}")
      @lin_equ[k+1][n+3] += m_f_k_xi*eta
      (-1..n+1).each do |i|
        f_i_xi = f_i(xi, i)
        @lin_equ[k+1][i+1] += m_f_k_xi * f_i_xi
        # puts("(#{i}, #{k}): f_#{i}(#{xi})=#{f_i_xi} p=#{m_f_k_xi * f_i_xi}")
      end
    end
    @count += m
    # puts "aggregate(xi=#{xi} eta=#{eta}) lin_equ=#{@lin_equ}"
  end

  def calc_cub
    @cub_sums = Array.new(n+3, 0)
    i = -1
    @cub = @lin_equ.map do |line|
      k = -1
      cline = line.map do |element|
        c = element.abs**3
        if (k <= n+1)
          @cub_sums[i+1]+=c
        end
        k += 1
        c
      end
      i += 1
      cline
    end
    # puts "lin_equ=#{@lin_equ}"
    # puts "cub_sums=#{@cub_sums}"
    # puts "cub=#{@cub}"
    STDOUT.flush
  end
  
  def solve
    calc_cub()
    (-1..n+1).each do |k|
      pivot_i = -1
      pivot_crit = 0
      (k..n+1).each do |i|
        if (@cub_sums[i+1] == 0)
          err = "cub_sums=#{@cub_sums} k=#{k} i=#{i} n=#{n} lin_equ=#{@lin_equ}"
          STDERR.puts "ERROR: #{err}"
          raise ZeroDivisionError.new(err)
        end
        pc = @cub[i+1][k+1] / @cub_sums[i+1]
        if (pc > pivot_crit)
          pivot_i = i
          pivot_crit = pc
        end
      end
      if (pivot_i != k)
        line = @lin_equ[pivot_i+1]
        @lin_equ[pivot_i+1] = @lin_equ[k+1]
        @lin_equ[k+1] = line
      end
      (k+1..n+1).each do |i|
        pivot = @lin_equ[k+1][k+1]
        if (pivot.kind_of? Integer)
          pivot = pivot.to_f
        end
        if (@lin_equ[i+1][k+1] != 0)
          q = @lin_equ[i+1][k+1]/pivot
          @lin_equ[i+1][k+1] = 0
          (k+1..n+2).each do |l|
            @lin_equ[i+1][l+1] -= q*@lin_equ[k+1][l+1]
          end
        end
      end
      # inefficient but ok as first step
      calc_cub
    end

    @results = Array.new(n+3, 0)
    (n+1).downto(-1) do |k|
      result = @lin_equ[k+1][n+3]
      # puts "k=#{k} result=#{result} (p)"
      (k+1..n+1).each do |l|
        result -= @results[l+1]*@lin_equ[k+1][l+1]
        # puts "k=#{k} l=#{l} result=#{result} (s)"
      end
      # puts "k=#{k} result=#{result} (q)"
      result /= @lin_equ[k+1][k+1]
      # puts "k=#{k} result=#{result}"
      @results[k+1] = result
    end
    calculate_coefficients()
    @results
  end

  def ser
    "coefficients=#{@coefficients}\ncount={@count}\nh={@h}\nn=#{@n}\nresults=#{@results}\nx_min=#{@x_min}\nx_max={@x_max}\nx=#{@x}\nx_i=#{@x_i}\n";
  end

  def self.deser(s)
    lines = s.split("\n")
    x_min = nil
    x_max = nil
    rmi=/x_min=(?<V>\S+)/
    rma=/x_max=(?<V>\S+)/
    lines.each do |line|
      m = line.match(rmi)
      if (m) then
        x_min = m[:V].to_f
      end
      m = line.match(rma)
      if (m) then
        x_max = m[:V].to_f
      end
    end
    if ((x_min.kind_of? Float) && (x_max.kind_of? Float)) then
      result = SplineApproximation.new(x_min, x_max)
      result.deser_internal(lines)
      result
    else
      raise "x_min or x_max missing"
    end
  end
    
  def deser_internal(lines)
    r_count = /count=(?<V>\d+)/
    r_h = /h=(?<V>\S+)/
    r_n = /n=(?<V>\d+)/
    r_results = /results=\[(?<V>.+)\]/
    r_x = /x=\[(?<V>.+)\]/
    r_x_i = /x_i=\[(?<V>.+)\]/
    lines.each do |line|
      m = line.match(r_count)
      if (m)
        @count = m[:V].to_i
      end
      m = line.match(r_h)
      if (m)
        @h = m[:V].to_f
      end
      m = line.match(r_n)
      if (m)
        @n = m[:V].to_i
      end
      m = line.match(r_results)
      if (m)
         @results = m[:V].split(/,/).map do |e|
          e.to_f
        end
      end
      m = line.match(r_x)
      if (m)
         @x = m[:V].split(/,/).map do |e|
          e.to_f
        end
      end
      m = line.match(r_x_i)
      if (m)
         @x_i = m[:V].split(/,/).map do |e|
          e.to_f
        end
      end
    end
    calculate_coefficients
  end
    
  def calculate_coefficients
    @coefficients = Array.new(n+6) { Array.new(4,0) }
    h=@h
    ext_result = ([0.0]*4).concat(@results).concat([0.0]*4)
    (-3..(n+2)).each do |i|
      a = ext_result[5+i-1]
      b = ext_result[5+i]
      c = ext_result[5+i+1]
      d = ext_result[5+i+2]
      @coefficients[i+3][0] = (-d+5*c+17*b+27*a)/8.0
      @coefficients[i+3][1] = (3*d-9*c+33*b-27*a)/(4.0*h)
      @coefficients[i+3][2] = (-3*d+15*c-21*b+9*a)/(2.0*h*h)
      @coefficients[i+3][3] = (d-3*c+3*b-a)/(1.0*h*h*h)
    end
  end

  # simple implementation of g
  def gs(x)
    (-1..n+1).map do |i|
      @results[i+1]*f_i(x, i)
    end.reduce(0) do |s, t|
      s+t
    end
  end

  # efficient implementation of g
  def g(x)
    if (x <= @x_minus_3)
      0
    elsif (x >= @x_n_plus_3)
      0
    else
      i = ((x - @x_min)/@h).floor
      xx = (x_i(i) + x_i(i+1))/2
      y = x-xx
      ((@coefficients[i+3][3]*y+@coefficients[i+3][2])*y+@coefficients[i+3][1])*y+@coefficients[i+3][0]
    end
  end
  
end
