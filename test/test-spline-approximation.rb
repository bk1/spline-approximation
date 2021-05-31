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

  def _test_aggregate_const_0
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

  def test_aggregate_const_1_m2
    check_aggregate_const(0, 1, 2, 1, 8)
    check_aggregate_const(0, 1, 2, 1, 7)
    check_aggregate_const(0, 1, 2, 1, 6)
    check_aggregate_const(0, 1, 2, 1, 5)
  end
  
  def test_aggregate_const_1_m10
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

  def test_aggregate_const_10_mx
    check_aggregate_const(0, 1, 2, 10, 5)
    check_aggregate_const(0, 1, 3, 10, 6)
    check_aggregate_const(0, 1, 4, 10, 7)
    check_aggregate_const(0, 1, 5, 10, 8)
    check_aggregate_const(0, 1, 6, 10, 9)
    check_aggregate_const(0, 1, 7, 10, 10)
    check_aggregate_const(0, 1, 90, 10, 93)
  end

  def check_aggregate_const(x_min, x_max, n, c, m)
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
      # puts("x=#{x} y=#{y} i=#{i} x_min=#{x_min} x_max=#{x_max} n=#{n} c=#{c} m=#{m}")
      # STDOUT.flush
      assert_in_delta(c, y, delta, "x=#{x} y=#{y} i=#{i} x_min=#{x_min} x_max=#{x_max} n=#{n} c=#{c} m=#{m}")
    end
    puts
  end
  
    
  end

# end of file testlongdecimal.rb
