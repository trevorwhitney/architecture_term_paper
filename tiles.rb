require 'rubygems'
require 'gruff'
require 'csv'

class TileResults
  attr_reader :results

  def initialize(file)
    @results = Array.new(5)
    5.times do |i|
      @results[i] = Array.new
    end

    CSV.foreach(file, :headers => :first_row, :return_headers => false) do |row|
      index = row[1].to_i/500 - 1
      @results[index] << row[3].to_f
    end
  end

  def graph(title, filename)
    g = Gruff::Line.new
    g.title = title

    g.title_font_size = 26
    g.legend_font_size = 16

    g.theme_37signals

    g.data("1500", @results[2])
    g.data("2000", @results[3])
    g.data("2500", @results[4])

    values = (16..64).to_a
    labels = {}

    values.each_with_index do |value, i|
      if i % 8 == 0
        labels[i] = value.to_s
      end
    end

    g.labels = labels
    g.hide_dots = true

    g.x_axis_label = "Tiles"
    g.y_axis_label = "Time (Seconds)"

    g.write('plots/' + filename)
  end

end

results = TileResults.new('tiles_final/tiles.csv')
puts results.results[4][0]
results.graph("Optimal Tile Size", 'tiles.png')