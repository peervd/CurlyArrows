from operate_analysis import analyze

'''

This script can be run to get custom feedback from your openAI account on a 
reaction mechanism. You can draw your reaction mechanism in ChemDoodle. To run
the analysis you have to copy your openai key, the exercise number and the
ChemDoodle drawing as JSON code. 

'''

# Copy the openAI key between quotation marks
# find your key under "specify later"
ai_key = False

# integer number which corresponds to the exercise number of the reaction mechanism
exercise = 1


# copy json_code from your reaction mechanism in ChemDoodle
json_code = '''{"m":[{"a":[{"x":45.6340373788567,"y":270.52272727272725,"i":"a0"},{"x":62.95454545454547,"y":260.52272727272725,"i":"a1"},{"x":80.27505353023423,"y":270.52272727272725,"i":"a2"},{"x":62.95454545454547,"y":240.52272727272728,"i":"a3","l":"O"},{"x":97.595561605923,"y":260.52272727272725,"i":"a4","l":"H"}],"b":[{"b":0,"e":1,"i":"b0"},{"b":1,"e":2,"i":"b1"},{"b":1,"e":3,"i":"b2","o":2},{"b":2,"e":4,"i":"b3"}]},{"a":[{"x":137.67949192431126,"y":260.29545454545456,"i":"a5","l":"O","c":-1}]},{"a":[{"x":241.0885828334022,"y":262.56818181818187,"i":"a6"},{"x":258.409090909091,"y":252.56818181818184,"i":"a7"},{"x":275.7295989847797,"y":262.56818181818187,"i":"a8"},{"x":258.409090909091,"y":232.56818181818184,"i":"a9","l":"O","c":-1}],"b":[{"b":0,"e":1,"i":"b4"},{"b":1,"e":2,"i":"b5","o":2},{"b":1,"e":3,"i":"b6"}]},{"a":[{"x":339.90443839407703,"y":260.1500000000001,"i":"a10"},{"x":357.22494646976577,"y":250.15000000000003,"i":"a11"},{"x":374.5454545454544,"y":260.1500000000001,"i":"a12"},{"x":357.22494646976577,"y":230.15000000000003,"i":"a13","l":"O"}],"b":[{"b":0,"e":1,"i":"b7"},{"b":1,"e":2,"i":"b8"},{"b":1,"e":3,"i":"b9","o":2}]},{"a":[{"x":472.906764651584,"y":262.56818181818187,"i":"a14"},{"x":490.22727272727275,"y":252.56818181818184,"i":"a15"},{"x":507.5477808029615,"y":262.56818181818187,"i":"a16"},{"x":490.22727272727275,"y":232.56818181818184,"i":"a17","l":"O"},{"x":524.8682888786502,"y":252.56818181818187,"i":"a18"},{"x":542.188796954339,"y":262.56818181818187,"i":"a19"},{"x":524.8682888786502,"y":232.56818181818187,"i":"a20","l":"O","c":-1},{"x":524.8682888786502,"y":272.56818181818187,"i":"a21"}],"b":[{"b":0,"e":1,"i":"b10"},{"b":1,"e":2,"i":"b11"},{"b":1,"e":3,"i":"b12","o":2},{"b":2,"e":4,"i":"b13"},{"b":4,"e":5,"i":"b14"},{"b":4,"e":6,"i":"b15"},{"b":4,"e":7,"i":"b16"}]},{"a":[{"x":577.4522191970385,"y":250.06818181818184,"i":"a22","l":"O"},{"x":594.7727272727273,"y":240.06818181818184,"i":"a23","l":"H"}],"b":[{"b":0,"e":1,"i":"b17"}]},{"a":[{"x":685.0839303183884,"y":258.2772727272728,"i":"a24"},{"x":702.4044383940771,"y":248.27727272727276,"i":"a25"},{"x":719.7249464697659,"y":258.2772727272728,"i":"a26"},{"x":737.0454545454546,"y":248.27727272727276,"i":"a27"},{"x":754.3659626211434,"y":258.2772727272728,"i":"a28"},{"x":737.0454545454546,"y":268.27727272727276,"i":"a29"},{"x":702.4044383940771,"y":228.27727272727276,"i":"a30","l":"O"},{"x":719.7249464697659,"y":238.27727272727276,"i":"a31","l":"O"}],"b":[{"b":0,"e":1,"i":"b18"},{"b":1,"e":2,"i":"b19"},{"b":2,"e":3,"i":"b20"},{"b":3,"e":4,"i":"b21"},{"b":3,"e":5,"i":"b22"},{"b":1,"e":6,"i":"b23"},{"b":1,"e":7,"i":"b24"},{"b":7,"e":3,"i":"b25"}]}],"s":[{"i":"s0","t":"Pusher","o1":"a5","o2":"a4","e":2},{"i":"s1","t":"Pusher","o1":"b3","o2":"b1","e":2},{"i":"s2","t":"Pusher","o1":"b2","o2":"a3","e":2},{"i":"s3","t":"Line","x1":172.906764651584,"y1":252.3409090909091,"x2":211.54312828794764,"y2":253.47727272727275,"a":"synthetic"},{"i":"s4","t":"Pusher","o1":"a9","o2":"b6","e":2},{"i":"s5","t":"Pusher","o1":"b5","o2":"a11","e":2},{"i":"s6","t":"Pusher","o1":"b9","o2":"a13","e":2},{"i":"s7","t":"Line","x1":406.99767374249313,"y1":244.38636363636365,"x2":449.05848175004246,"y2":244.38636363636365,"a":"synthetic"},{"i":"s8","t":"Line","x1":617.7272796630859,"y1":246.5909090909091,"x2":645.0236709118684,"y2":246.5909090909091,"a":"synthetic"},{"i":"s9","t":"Pusher","o1":"a20","o2":"a15","e":2},{"i":"s10","t":"Pusher","o1":"a15","o2":"a17","e":2},{"i":"s11","t":"Pusher","o1":"a17","o2":"a23","e":2},{"i":"s12","t":"Pusher","o1":"b17","o2":"a22","e":2}]}'''

''' now run the script '''
output = analyze(ai_key,exercise,json_code)

print(output)