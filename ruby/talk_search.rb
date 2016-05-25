require 'net/http'
require 'json'

# A Ruby script that returns information from Talk via the Zooniverse Ouroboros API
#
# Originally written by Michael Parrish; adapted by Kyle Willett

def search(project: nil, text: nil, kind: 'subject', tags: { }, page: 1, results: [])
  raise 'no project specified' unless project
  tag_query = tags.each_pair.collect{ |k, v| "tags[#{ k }]=#{ v }" }.join '&'
  uri = URI.parse "https://api.zooniverse.org/projects/#{ project }/talk/search?text=#{ text }&kind=#{ kind }&#{ tag_query }&per_page=20&page=#{ page }"
  req = Net::HTTP::Get.new uri.to_s
  http = Net::HTTP.new uri.host, uri.port
  http.use_ssl = true
  res = http.request req
  json = JSON.parse res.body
  
  pages = (json['total'] / json['per_page'].to_f).ceil
  
  # More than 1,000 results
  if page == 1 && pages > 50
    puts "\n\nThis query has #{ json['total'] } results."
    puts "It could take a long time and degrade server performance."
    puts "Are you really sure you want to run this query? (y/n)"
    return unless gets.chomp =~ /y/i
  end
  
  if json['page'] < pages
    puts "#{ json['page'] } / #{ pages }"
    search project: project, text: text, kind: kind, tags: tags, page: page + 1, results: results + json['results']
  else
    results + json['results']
  end
end

=begin

# Below are some example searches

# search for subjects tagged with both merger and spiral
merger_and_spiral_subjects = search project: 'galaxy_zoo', tags: { merger: true, spiral: true }

# search for subjects tagged with both merger and spiral, but not overlap and not disturbed (all conditions are ANDed together)
mixed_booleans = search project: 'galaxy_zoo', tags: { merger: true, spiral: true, overlap: false, disturbed: false }

# search for other types, like 'collection', or 'discussion'
merger_collections = search project: 'galaxy_zoo', kind: 'collection', tags: { merger: true, spiral: true }

# text search is combinable as well
pretty_mergers = search project: 'galaxy_zoo', text: 'pretty', tags: { merger: true }
=end

# Examples: queries for giant radio galaxies. Post-processed in Python (talk_search.py)

radio_giant = search project: 'radio', text: 'giant'
radio_large = search project: 'radio', text: 'large'
radio_huge = search project: 'radio', text: 'huge'
radio_kpc = search project: 'radio', text: 'kpc'
radio_mpc = search project: 'radio', text: 'Mpc'
radio_overedge = search project: 'radio', text: 'overedge'

# write the output to a file if you'd rather work with in another language
File.open('talk_searches/radio_giant.json', 'w'){ |out| out.puts JSON.dump(radio_giant) }
File.open('talk_searches/radio_large.json', 'w'){ |out| out.puts JSON.dump(radio_large) }
File.open('talk_searches/radio_huge.json', 'w'){ |out| out.puts JSON.dump(radio_huge) }
File.open('talk_searches/radio_kpc.json', 'w'){ |out| out.puts JSON.dump(radio_kpc) }
File.open('talk_searches/radio_mpc.json', 'w'){ |out| out.puts JSON.dump(radio_mpc) }
File.open('talk_searches/radio_overedge.json', 'w'){ |out| out.puts JSON.dump(radio_overedge) }

# example response structure
example_result = <<-JSON
[
  {
    "id": "AGZ0003wd2",
    "tags": [
      "merger",
      "tidal",
      "dustlanes",
      "dustlane",
      "tidaltails",
      "mrrger",
      "tidaldebris",
      "verynice",
      "collision",
      "notam_wrong",
      "polarring"
    ],
    "updated_at": "2013-03-06T23:03:25Z",
    "location": {
      "thumbnail": "http://www.galaxyzoo.org.s3.amazonaws.com/subjects/thumbnail/1237672026246545455.jpg",
      "standard": "http://www.galaxyzoo.org.s3.amazonaws.com/subjects/standard/1237672026246545455.jpg",
      "inverted": "http://www.galaxyzoo.org.s3.amazonaws.com/subjects/inverted/1237672026246545455.jpg"
    },
    "name": "AGZ0003wd2",
    "project_id": "502a90cd516bcb060c000001",
    "kind": "subject",
    "_score": 69.35612,
    "_type": "focus",
    "_index": "talk",
    "_version": null,
    "sort": null,
    "highlight": null,
    "_explanation": null
  }
]
JSON
