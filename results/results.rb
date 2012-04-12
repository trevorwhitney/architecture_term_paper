require 'rubygems'
require 'gruff'
require 'pry'

class Results

  def initialize( o0_results, 
                  o1_results, 
                  o2_results, 
                  o3_results, 
                  o3_funroll_results)

    @o0_results = o0_results
    @o1_results = o1_results
    @o2_results = o2_results
    @o3_results = o3_results
    @o3_funroll_results = o3_funroll_results

    @best_times = [ { :algorithm => "", :time => 200},
                  { :algorithm => "", :time => 200},
                  { :algorithm => "", :time => 200},
                  { :algorithm => "", :time => 200},
                  { :algorithm => "", :time => 200} ]
  end

  def best_time(results, algorithm)
    @best_times.each_with_index do |best, i|
      if results.results[i+2][9] < best[:time]
        @best_times[i][:algorithm] = algorithm
        @best_times[i][:time] = results.results[i+2][9]
      end
    end
  end


  def generate_results
    results = ResultsFile.new(@o0_results)
    best_time(results, "o0")
    results.graph("No Compiler Optimizations", "O0.png")

    results = ResultsFile.new(@o1_results)
    best_time(results, "o1")
    results.graph("O1 Compiler Optimization", "O1.png")

    results = ResultsFile.new(@o2_results)
    best_time(results, "o2")
    results.graph("O2 Compiler Optimization", "O2.png")

    results = ResultsFile.new(@o3_results)
    best_time(results, "o3")
    results.graph("O3 Compiler Optimization", "O3.png")

    results = ResultsFile.new(@o3_funroll_results)
    best_time(results, "o3 w/ funroll")
    results.graph("O3 & Funroll Compiler Optimizations", "O3_funroll.png")

    graph_overall_performance
  end

  def graph_overall_performance
    o1_best = @best_times[0][:time]*1.0
    o2_best = @best_times[1][:time]*1.0
    o3_best = @best_times[2][:time]*1.0
    o4_best = @best_times[3][:time]*1.0
    o5_best = @best_times[4][:time]*1.0

    gflop = (2*5000**3)/10**9

    o1_gflop = gflop/o1_best
    o2_gflop = gflop/o2_best
    o3_gflop = gflop/o3_best
    o4_gflop = gflop/o4_best
    o5_gflop = gflop/o5_best

    o1_performance = (o1_gflop/13.6) * 100
    o2_performance = (o2_gflop/13.6) * 100
    o3_performance = (o3_gflop/13.6) * 100
    o4_performance = (o4_gflop/13.6) * 100
    o5_performance = (o5_gflop/13.6) * 100

    times = Array.new()
    labels = Array.new()

    algorithms = %w(step01 step02 step03 step04 step05)

    @best_times.each_with_index do |best, i|
      times << best[:time]
      labels << "#{algorithms[i]} #{best[:algorithm]}"
    end

    performance = [o1_performance, o2_performance, o3_performance, 
      o4_performance, o5_performance]

    gflops = [o1_gflop, o2_gflop, o3_gflop, o4_gflop, o5_gflop]

    write_performance("Best Execution Times", labels, times, 
      "Algorithm & Optimizations", "Time (Seconds)", 'plots/best_times.png', 
      50, 100)

    write_performance("Percent of Max Performance", labels, performance, 
      "Algorithm & Optimizations", "Percent of Max Performance", 
      'plots/percent_performance.png', 0, 50)

    write_performance("Performance in Gflops", labels, gflops, 
      "Algorithm & Optimizations", "Glops", 'plots/gflops.png', 0, 5)
  end


  def write_performance(title, labels, values, x, y, filename, min = nil, max = nil)
    g = Gruff::Bar.new
    g.title = title

    g.title_font_size = 26
    g.legend_font_size = 16

    g.theme_37signals

    labels.each_with_index do |label, i|
      g.data(label, values[i])
    end

    g.minimum_value = min unless min.nil?
    g.maximum_value = max unless max.nil?

    g.x_axis_label = x
    g.y_axis_label = y

    g.write(filename)
  end

end

