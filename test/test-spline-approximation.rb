#!/usr/bin/env ruby
#
# testlongdecimal.rb -- runit test for long-decimal.rb
#
# (C) Karl Brodowsky (IT Sky Consulting GmbH) 2006-2016
#
# TAG:       $TAG v1.00.04$
# Author:    $Author: bk1 $ (Karl Brodowsky)
#

require 'complex'

$test_type = nil
if ((RUBY_VERSION.match /^1\./) || (RUBY_VERSION.match /^2\.[01]/)) then
  require 'test/unit'
  $test_type = :v20
else
  require 'minitest/autorun'
  require 'test/unit/assertions'
  include Test::Unit::Assertions
  $test_type = :v21
end

load "lib/spline-approximation.rb"

if ($test_type == :v20)
  class UnitTest < Test::Unit::TestCase
  end
else
  class UnitTest < MiniTest::Test
  end
end

#
# test class for SplineApproximation 
#
class TestSplineApproximation_class < UnitTest

  #
  # test split_to_words and merge_from_words
  #
  def test_initialize_good_params
    print "\ntest_initialize_good_params [#{Time.now}]: "
    check_initialize_good_params(0, 1, 1)
    check_initialize_good_params(0, 1, 2)
    check_initialize_good_params(0, 1, 10)
    check_initialize_good_params(0.0, 1.0, 1)
    check_initialize_good_params(0.0, 1.0, 2)
    check_initialize_good_params(0.0, 1.0, 10)
  end

  def check_initialize_good_params(x_min, x_max, n)
    s = SplineApproximation.new(x_min, x_max, n)
    assert_equal(x_min, s.x_min)
    assert_equal(x_max, s.x_max)
    assert_equal(n, s.n)
  end

  #
  # test split_to_words and merge_from_words
  #
  def test_initialize_bad_params
    print "\ntest_initialize_bad_params [#{Time.now}]: "
    check_initialize_bad_params(0, 0, 1)
    check_initialize_bad_params(0, 1, 0)
    check_initialize_bad_params(0, -1, 10)
    check_initialize_bad_params(Complex(1, 0), 1.0, 1)
    check_initialize_bad_params(0.0, Complex(1.0), 2)
    check_initialize_bad_params(Complex(0.0), Complex(1.0), 10)
  end

  def check_initialize_bad_params(x_min, x_max, n)
    e = assert_raises RuntimeError do
      s = SplineApproximation.new(x_min, x_max, n)
    end
    puts e
  end

  def test_aggregate_const_0
    check_aggregate_const(0, 1, 1, 0, 6)
    check_aggregate_const(0, 1, 1, 0, 5)
    check_aggregate_const(0, 1, 1, 0, 4)

    check_aggregate_const(0, 1, 2, 0, 7)
    check_aggregate_const(0, 1, 2, 0, 6)
    check_aggregate_const(0, 1, 2, 0, 5)

    check_aggregate_const(0, 1, 10, 0, 15)
    check_aggregate_const(0, 1, 10, 0, 14)
    check_aggregate_const(0, 1, 10, 0, 13)
  end

  def test_aggregate_const_1
    check_aggregate_const(0, 1, 1, 1, 8)
    check_aggregate_const(0, 1, 1, 1, 7)
    check_aggregate_const(0, 1, 1, 1, 6)
    check_aggregate_const(0, 1, 1, 1, 5)
    check_aggregate_const(0, 1, 1, 1, 4)
  end

  def test_aggregate_const_1_n2
    check_aggregate_const(0, 1, 2, 1, 8)
    check_aggregate_const(0, 1, 2, 1, 7)
    check_aggregate_const(0, 1, 2, 1, 6)
    check_aggregate_const(0, 1, 2, 1, 5)
  end
  
  def test_aggregate_const_1_n10
    check_aggregate_const(0, 1, 10, 1, 16)
    check_aggregate_const(0, 1, 10, 1, 15)
    check_aggregate_const(0, 1, 10, 1, 14)
    check_aggregate_const(0, 1, 10, 1, 13)
  end

  def test_aggregate_const_10
    check_aggregate_const(0, 1, 1, 10, 8)
    check_aggregate_const(0, 1, 1, 10, 7)
    check_aggregate_const(0, 1, 1, 10, 6)
    check_aggregate_const(0, 1, 1, 10, 5)
    check_aggregate_const(0, 1, 1, 10, 4)
  end

  def test_aggregate_const_10_nx
    check_aggregate_const(0, 1, 2, 10, 5)
    check_aggregate_const(0, 1, 3, 10, 6)
    check_aggregate_const(0, 1, 4, 10, 7)
    check_aggregate_const(0, 1, 5, 10, 8)
    check_aggregate_const(0, 1, 6, 10, 9)
    check_aggregate_const(0, 1, 7, 10, 10)
    check_aggregate_const(0, 1, 90, 10, 93)
  end

  def check_aggregate_const(x_min, x_max, n, c, m)
    puts "------------------------------------------------------------"
    delta = c.to_f / 3
    s = SplineApproximation.new(x_min, x_max, n)
    (0..m).each do |i|
      xi = (x_min * (m-i) + x_max * i).to_f / m
      s.aggregate(xi, c, 1)
    end
    s.solve
    puts "x_min=#{x_min} x_max=#{x_max} n=#{n} c=#{c} m=#{m} r=#{s.results}"
    STDOUT.flush
    (0..2*m).each do |i|
      x = (x_min * (2*m-i) + x_max * i).to_f / (2*m)
      y = s.g(x)
      ys = s.gs(x)
      puts("x=#{x} y=#{y} i=#{i} x_min=#{x_min} x_max=#{x_max} n=#{n} c=#{c} m=#{m}")
      assert_in_delta(c, y, delta, "x=#{x} y=#{y} i=#{i} x_min=#{x_min} x_max=#{x_max} n=#{n} c=#{c} m=#{m}")
      assert_in_delta(c, ys, delta, "x=#{x} y=#{y} ys=#{ys} i=#{i} x_min=#{x_min} x_max=#{x_max} n=#{n} c=#{c} m=#{m}")
    end
    puts "------------------------------------------------------------"
    puts
  end


  def test_aggregate_linear_0
    check_aggregate_linear(0, 1, 1, 0, 6)
    check_aggregate_linear(0, 1, 1, 0, 5)
    check_aggregate_linear(0, 1, 1, 0, 4)

    check_aggregate_linear(0, 1, 2, 0, 7)
    check_aggregate_linear(0, 1, 2, 0, 6)
    check_aggregate_linear(0, 1, 2, 0, 5)

    check_aggregate_linear(0, 1, 10, 0, 15)
    check_aggregate_linear(0, 1, 10, 0, 14)
    check_aggregate_linear(0, 1, 10, 0, 13)
  end

  def test_aggregate_linear_1
    check_aggregate_linear(0, 1, 1, 1, 8)
    check_aggregate_linear(0, 1, 1, 1, 7)
    check_aggregate_linear(0, 1, 1, 1, 6)
    check_aggregate_linear(0, 1, 1, 1, 5)
    check_aggregate_linear(0, 1, 1, 1, 4)
  end

  def test_aggregate_linear_1_n2
    check_aggregate_linear(0, 1, 2, 1, 8)
    check_aggregate_linear(0, 1, 2, 1, 7)
    check_aggregate_linear(0, 1, 2, 1, 6)
    check_aggregate_linear(0, 1, 2, 1, 5)
  end
  
  def test_aggregate_linear_1_n10
    check_aggregate_linear(0, 1, 10, 1, 16)
    check_aggregate_linear(0, 1, 10, 1, 15)
    check_aggregate_linear(0, 1, 10, 1, 14)
    check_aggregate_linear(0, 1, 10, 1, 13)
  end

  def test_aggregate_linear_10
    check_aggregate_linear(0, 1, 1, 10, 8)
    check_aggregate_linear(0, 1, 1, 10, 7)
    check_aggregate_linear(0, 1, 1, 10, 6)
    check_aggregate_linear(0, 1, 1, 10, 5)
    check_aggregate_linear(0, 1, 1, 10, 4)
  end

  def test_aggregate_linear_10_nx
    check_aggregate_linear(0, 1, 2, 10, 5)
    check_aggregate_linear(0, 1, 3, 10, 6)
    check_aggregate_linear(0, 1, 4, 10, 7)
    check_aggregate_linear(0, 1, 5, 10, 8)
    check_aggregate_linear(0, 1, 6, 10, 9)
    check_aggregate_linear(0, 1, 7, 10, 10)
    check_aggregate_linear(0, 1, 90, 10, 93)
  end

  def check_aggregate_linear(x_min, x_max, n, y0, m, step = 1)
    puts "------------------------------------------------------------"
    delta = (y0.abs + (y0+m*step).abs).to_f / 1000
    s = SplineApproximation.new(x_min, x_max, n)
    (0..m).each do |i|
      xi = (x_min * (m-i) + x_max * i).to_f / m
      s.aggregate(xi, y0+i*step, 1)
    end
    s.solve
    # c = s.coefficients
    # c.each do |l|
    #   p(l)
    # end
    puts "x_min=#{x_min} x_max=#{x_max} n=#{n} y0=#{y0} step=#{step} m=#{m} r=#{s.results}"
    STDOUT.flush
    (0..2*m).each do |i|
      x = (x_min * (2*m-i) + x_max * i).to_f / (2*m)
      y = s.g(x)
      ys = s.gs(x)
      ye = y0 + 0.5*i*step
      puts("x=#{x} y=#{y} ys=#{ys} ye=#{ye} i=#{i} x_min=#{x_min} x_max=#{x_max} n=#{n} ye=#{ye} m=#{m}")
      # assert_in_delta(ye, y, delta, "x=#{x} y=#{y} ys=#{ys} i=#{i} x_min=#{x_min} x_max=#{x_max} n=#{n} ye=#{ye} m=#{m}")
      assert_in_delta(ye, ys, delta, "x=#{x} y=#{y} ys=#{ys} i=#{i} x_min=#{x_min} x_max=#{x_max} n=#{n} ye=#{ye} m=#{m}")
    end
    puts "------------------------------------------------------------"
    puts
  end

  def test_aggregate_quadratic_0
    check_aggregate_quadratic(0, 1, 1, 0, 6)
    check_aggregate_quadratic(0, 1, 1, 0, 5)
    check_aggregate_quadratic(0, 1, 1, 0, 4)

    check_aggregate_quadratic(0, 1, 2, 0, 7)
    check_aggregate_quadratic(0, 1, 2, 0, 6)
    check_aggregate_quadratic(0, 1, 2, 0, 5)

    check_aggregate_quadratic(0, 1, 10, 0, 15)
    check_aggregate_quadratic(0, 1, 10, 0, 14)
    check_aggregate_quadratic(0, 1, 10, 0, 13)
  end

  def test_aggregate_quadratic_1
    check_aggregate_quadratic(0, 1, 1, 1, 8)
    check_aggregate_quadratic(0, 1, 1, 1, 7)
    check_aggregate_quadratic(0, 1, 1, 1, 6)
    check_aggregate_quadratic(0, 1, 1, 1, 5)
    check_aggregate_quadratic(0, 1, 1, 1, 4)
  end

  def test_aggregate_quadratic_1_n2
    check_aggregate_quadratic(0, 1, 2, 1, 8)
    check_aggregate_quadratic(0, 1, 2, 1, 7)
    check_aggregate_quadratic(0, 1, 2, 1, 6)
    check_aggregate_quadratic(0, 1, 2, 1, 5)
  end
  
  def test_aggregate_quadratic_1_n10
    check_aggregate_quadratic(0, 1, 10, 1, 16)
    check_aggregate_quadratic(0, 1, 10, 1, 15)
    check_aggregate_quadratic(0, 1, 10, 1, 14)
    check_aggregate_quadratic(0, 1, 10, 1, 13)
  end

  def test_aggregate_quadratic_10
    check_aggregate_quadratic(0, 1, 1, 10, 8)
    check_aggregate_quadratic(0, 1, 1, 10, 7)
    check_aggregate_quadratic(0, 1, 1, 10, 6)
    check_aggregate_quadratic(0, 1, 1, 10, 5)
    check_aggregate_quadratic(0, 1, 1, 10, 4)
  end

  def test_aggregate_quadratic_10_nx
    check_aggregate_quadratic(0, 1, 2, 10, 5)
    check_aggregate_quadratic(0, 1, 3, 10, 6)
    check_aggregate_quadratic(0, 1, 4, 10, 7)
    check_aggregate_quadratic(0, 1, 5, 10, 8)
    check_aggregate_quadratic(0, 1, 6, 10, 9)
    check_aggregate_quadratic(0, 1, 7, 10, 10)
    check_aggregate_quadratic(0, 1, 90, 10, 93)
  end

  def check_aggregate_quadratic(x_min, x_max, n, y0, m, step = 1)
    puts "------------------------------------------------------------"
    delta = (y0.abs + (y0+m*step).abs).to_f / 1000
    s = SplineApproximation.new(x_min, x_max, n)
    (0..m).each do |i|
      xi = (x_min * (m-i) + x_max * i).to_f / m
      qi = (y0+i*step) ** 2
      s.aggregate(xi, qi, 1)
    end
    s.solve
    # c = s.coefficients
    # c.each do |l|
    #   p(l)
    # end
    puts "x_min=#{x_min} x_max=#{x_max} n=#{n} y0=#{y0} step=#{step} m=#{m} r=#{s.results}"
    STDOUT.flush
    (0..2*m).each do |i|
      x = (x_min * (2*m-i) + x_max * i).to_f / (2*m)
      q = s.g(x)
      qs = s.gs(x)
      qe = (y0 + 0.5*i*step)**2
      puts("x=#{x} q=#{q} qs=#{qs} qe=#{qe} i=#{i} x_min=#{x_min} x_max=#{x_max} n=#{n} qe=#{qe} m=#{m}")
      # assert_in_delta(qe, q, delta, "x=#{x} q=#{q} qs=#{qs} i=#{i} x_min=#{x_min} x_max=#{x_max} n=#{n} qe=#{qe} m=#{m}")
      assert_in_delta(qe, qs, delta, "x=#{x} q=#{q} qs=#{qs} i=#{i} x_min=#{x_min} x_max=#{x_max} n=#{n} qe=#{qe} m=#{m}")
    end
    puts "------------------------------------------------------------"
    puts
  end

  def test_aggregate_cubic_0
    check_aggregate_cubic(0, 1, 1, 0, 6)
    check_aggregate_cubic(0, 1, 1, 0, 5)
    check_aggregate_cubic(0, 1, 1, 0, 4)

    check_aggregate_cubic(0, 1, 2, 0, 7)
    check_aggregate_cubic(0, 1, 2, 0, 6)
    check_aggregate_cubic(0, 1, 2, 0, 5)

    check_aggregate_cubic(0, 1, 10, 0, 15)
    check_aggregate_cubic(0, 1, 10, 0, 14)
    check_aggregate_cubic(0, 1, 10, 0, 13)
  end

  def test_aggregate_cubic_1
    check_aggregate_cubic(0, 1, 1, 1, 8)
    check_aggregate_cubic(0, 1, 1, 1, 7)
    check_aggregate_cubic(0, 1, 1, 1, 6)
    check_aggregate_cubic(0, 1, 1, 1, 5)
    check_aggregate_cubic(0, 1, 1, 1, 4)
  end

  def test_aggregate_cubic_1_n2
    check_aggregate_cubic(0, 1, 2, 1, 8)
    check_aggregate_cubic(0, 1, 2, 1, 7)
    check_aggregate_cubic(0, 1, 2, 1, 6)
    check_aggregate_cubic(0, 1, 2, 1, 5)
  end
  
  def test_aggregate_cubic_1_n10
    check_aggregate_cubic(0, 1, 10, 1, 16)
    check_aggregate_cubic(0, 1, 10, 1, 15)
    check_aggregate_cubic(0, 1, 10, 1, 14)
    check_aggregate_cubic(0, 1, 10, 1, 13)
  end

  def test_aggregate_cubic_10
    check_aggregate_cubic(0, 1, 1, 10, 8)
    check_aggregate_cubic(0, 1, 1, 10, 7)
    check_aggregate_cubic(0, 1, 1, 10, 6)
    check_aggregate_cubic(0, 1, 1, 10, 5)
    check_aggregate_cubic(0, 1, 1, 10, 4)
  end

  def test_aggregate_cubic_10_nx
    check_aggregate_cubic(0, 1, 2, 10, 5)
    check_aggregate_cubic(0, 1, 3, 10, 6)
    check_aggregate_cubic(0, 1, 4, 10, 7)
    check_aggregate_cubic(0, 1, 5, 10, 8)
    check_aggregate_cubic(0, 1, 6, 10, 9)
    check_aggregate_cubic(0, 1, 7, 10, 10)
    check_aggregate_cubic(0, 1, 90, 10, 93)
  end

  def check_aggregate_cubic(x_min, x_max, n, y0, m, step = 1)
    puts "------------------------------------------------------------"
    delta = (y0.abs + (y0+m*step).abs).to_f / 1000
    s = SplineApproximation.new(x_min, x_max, n)
    (0..m).each do |i|
      xi = (x_min * (m-i) + x_max * i).to_f / m
      ci = (y0+i*step) ** 2
      s.aggregate(xi, ci, 1)
    end
    s.solve
    # c = s.coefficients
    # c.each do |l|
    #   p(l)
    # end
    puts "x_min=#{x_min} x_max=#{x_max} n=#{n} y0=#{y0} step=#{step} m=#{m} r=#{s.results}"
    STDOUT.flush
    (0..2*m).each do |i|
      x = (x_min * (2*m-i) + x_max * i).to_f / (2*m)
      c = s.g(x)
      cs = s.gs(x)
      ce = (y0 + 0.5*i*step)**2
      puts("x=#{x} c=#{c} cs=#{cs} ce=#{ce} i=#{i} x_min=#{x_min} x_max=#{x_max} n=#{n} ce=#{ce} m=#{m}")
      # assert_in_delta(ce, c, delta, "x=#{x} c=#{c} cs=#{cs} i=#{i} x_min=#{x_min} x_max=#{x_max} n=#{n} ce=#{ce} m=#{m}")
      assert_in_delta(ce, cs, delta, "x=#{x} c=#{c} cs=#{cs} i=#{i} x_min=#{x_min} x_max=#{x_max} n=#{n} ce=#{ce} m=#{m}")
    end
    puts "------------------------------------------------------------"
    puts
  end
    
end

# end of file testlongdecimal.rb
