
report_name=$1

sha=$(jq -r '.sha' $report_name)
post_adress=https://api.github.com/repos/qjutard/chl_bbp_ttt/statuses/$sha

echo $post_adress

finished=$(jq '{"state":.finished,"description":.finished_description,"context":"finished"}' $report_name)
no_errors=$(jq '{"state":.no_errors,"description":.no_errors_description,"context":"no_errors"}' $report_name)
filechecker=$(jq '{"state":.filechecker,"description":.filechecker_description,"context":"filechecker"}' $report_name)

echo $finished
echo $no_errors
echo $filechecker

curl -u qjutard -X POST -H 'Content-Type: application/json' --data "$finished" $post_adress
curl -u qjutard -X POST -H 'Content-Type: application/json' --data "$no_errors" $post_adress
curl -u qjutard -X POST -H 'Content-Type: application/json' --data "$filechecker" $post_adress
