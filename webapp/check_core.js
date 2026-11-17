const fs=require('fs'); const M=require('./core.js');
const gs=JSON.parse(fs.readFileSync('genesets.json')).genesets; const ck=JSON.parse(fs.readFileSync('clock.json'));
const c=M.parseCounts(fs.readFileSync('iss_DA_counts.csv','utf8'));
const {scores}=M.scoreAll(c,gs,ck);
const exp=M.parseCounts(fs.readFileSync('iss_expected_scores.csv','utf8'));
let maxErr=0; exp.genes.forEach((s,i)=>{const j=c.samples.indexOf(s); exp.samples.forEach((p,k)=>{maxErr=Math.max(maxErr,Math.abs(scores[p][j]-exp.X[i][k]));});});
console.log('programs compared:',exp.samples.join(','));
console.log('JS vs Python max abs diff (real ISS data):', maxErr);