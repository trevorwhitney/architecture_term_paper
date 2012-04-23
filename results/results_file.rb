require 'rubygems'
require 'gruff'
require 'csv'

class ResultsFile
  attr_reader :results
  attr_reader :error

  def initialize(file)
    @results = Array.new
    @error = Array.new
    7.times do
      @results << Array.new(10)
      @error << Array.new(10)
    end
    o0 = CSV.foreach(file, :header => :first_row, :return_headers => false) do |row|
      index = row[1].to_i/500 - 1
      
      if row[0] == "atlas"
        parse_data(0, index, row)
      elsif row[0] == "my_dgemm"
        parse_data(1, index, row)
      elsif row[0] == "step01"
        parse_data(2, index, row)
      elsif row[0] == "step02"
        parse_data(3, index, row)
      elsif row[0] == "step03"
        parse_data(4, index, row)
      elsif row[0] == "step04"
        parse_data(5, index, row)
      elsif row[0] == "step05"
        parse_data(6, index, row)
      end

    end
  end

  def parse_data(algorithm_id, index, row)
    @results[algorithm_id][index] = row[3].to_f
    @error[algorithm_id][index] = row[4].to_f
  end

  def graph(title, filename)
    no_dgemm_filename = File.basename(filename, '.png') + "_noDgemm.png"
    no_dgemm_title = title + " (w/o my_dgemm)"
    
    error_filename = File.basename(filename, '.png') + "_error.png"
    error_title = title + " Error"

    results_graph(title, filename)
    results_graph(no_dgemm_title, no_dgemm_filename, false)

    error_graph(error_title, error_filename)
    #puts @error[6]
  end

  def results_graph(title, filename, dgemm = true)
    g = Gruff::Line.new
    g.title = title

    g.title_font_size = 26
    g.legend_font_size = 16

    g.theme_37signals

    g.data("Atlas", @results[0])
    g.data("Step01", @results[2])
    g.data("Step02", @results[3])
    g.data("Step03", @results[4])
    g.data("Step04", @results[5])
    g.data("Step05", @results[6])
    g.data("My Dgemm", @results[1]) if dgemm

    g.labels = {0 => '500', 2 => '1500', 4 => '2500', 6 => '3500', 8 => '4500',
      9 => '5000'}

    g.x_axis_label = "Size of Matrix (N)"
    g.y_axis_label = "Time (Seconds)"

    g.write('plots/' + filename)
  end

  def error_graph(title, filename)
    g = Gruff::Bar.new
    g.title = title

    g.title_font_size = 26
    g.legend_font_size = 16

    g.theme_37signals

    g.data("Atlas", @error[0])
    g.data("My Dgemm", @error[1])
    g.data("Step01", @error[2])
    g.data("Step02", @error[3])
    g.data("Step03", @error[4])
    g.data("Step04", @error[5])
    g.data("Step05", @error[6])

    g.labels = {0 => '500', 2 => '1500', 4 => '2500', 6 => '3500', 8 => '4500',
      9 => '5000'}

    g.x_axis_label = "Size of Matrix (N)"
    g.y_axis_label = "Error"

    g.write('plots/' + filename)
  end

end


