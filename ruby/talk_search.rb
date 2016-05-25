require 'net/http'
require 'json'

# A Ruby script that returns information from Talk via the Zooniverse Ouroboros API. 
# Results are saved as JSON files.
#
# Originally written by Michael Parrish (Adler/Zooniverse); adapted by Kyle Willett (UMN)
#
# Example:
#
# >> ruby talk_search.rb

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

########################################
# Search parameters
#
# This is the query that's actually executed. Change the project and/or text if you want to do a different search.
#
########################################

radio_giant = search project: 'radio', text: 'giant'

# Write the output to a file if you'd rather analyze the results in another language (like Python)
# Post-processed in Python (talk_search.py)

File.open('talk_searches/radio_giant.json', 'w'){ |out| out.puts JSON.dump(radio_giant) }

# Example response structure
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
