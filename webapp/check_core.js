const fs=require('fs'); const M=require('./core.js');
const gs=JSON.parse(fs.readFileSync('genesets.json')).genesets;
const c=M.parseCounts(fs.readFileSync('demo_counts.csv','utf8'));
const {scores}=M.scoreAll(c,gs);
const exp=M.parseCounts(fs.readFileSync('demo_expected_scores.csv','utf8')); // samples x programs, read as table
let maxErr=0; const progs=exp.samples;
exp.genes.forEach((s,i)=>{const j=c.samples.indexOf(s); progs.forEach((p,k)=>{maxErr=Math.max(maxErr,Math.abs(scores[p][j]-exp.X[i][k]));});});
console.log('JS vs Python max abs diff:', maxErr);
