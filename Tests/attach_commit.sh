
test_log=$1

sha=$(git rev-parse HEAD)

git log -1

read -p "
This commit will be associated to test log $1, confirm (y/n) ?" choice
case "$choice" in
	y|Y ) echo "continuing...";;
	n|N ) exit 0;;
	*) echo "invalid choice" && exit 1;;
esac
 
new_report=$(jq --arg sha "$sha" '.sha=$sha' $1)

echo $new_report | jq '.' > $1


