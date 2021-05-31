# encoding: utf-8

#
# Rakefile for spline-approximation project
#
# (C) Karl Brodowsky (IT Sky Consulting GmbH) 2021
#
# Author:    bk1 (Karl Brodowsky)
#

# require "bundler/gem_tasks"
require 'rubygems'
require 'bundler'

require 'rake/testtask'
Rake::TestTask.new(:test) do |test|
  test.libs << 'lib' << 'test'
  test.pattern = 'test/test*.rb'
  test.verbose = true
end

require 'rdoc/task'
Rake::RDocTask.new do |rdoc|
  version = '0.1' # File.exist?('VERSION') ? File.read('VERSION').strip : LongDecimalSupport::VERSION

  rdoc.rdoc_dir = 'rdoc'
  rdoc.title = "spline-approximation #{version}"
  rdoc.main = "README"
  rdoc.rdoc_files.include('README*')
  rdoc.rdoc_files.include('lib/*.rb')
end

# end of file Rakefile

