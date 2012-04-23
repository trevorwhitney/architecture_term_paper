values = (16..64).to_a

labels = {}

values.each_with_index do |value, i|
  labels[i] = value.to_s
end

puts labels