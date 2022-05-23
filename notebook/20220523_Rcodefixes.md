# Fixes needed to the R code

* `final.df$cor_af_sal_noutliers` and `temp`: it would be better to output the number of QTN outliers and the number of truly neutral outliers - I think this calculation migth be wrong

* `LEA3.2_lfmm2_num_causal_sig_temp` there is a second variable `LEA3.2_lfmm2_num_causal_sig_temp` and analogous variable for salinity seems to be missing `LEA3.2_lfmm2_num_neut_sig_sal`


<head><style>
        table {
              font-family: times ;
color:  black ;
text-align: right;}
        th {
              padding: 1px 1px 5px 5px;
	        }
        td {
             padding: 1px 1px 5px 5px; }
      </style></head><table align="center" style="border-collapse: collapse; caption-side:top; font-size:11pt;"><caption style="text-align:center;"></caption><tr>
<th style="border-left: 0px solid black;background-color: #FFFFFF;border-top: 2px solid gray;border-bottom: 1px solid gray;">&nbsp;</th>
<th <th align="center" style="font-weight: normal;border-left: 0px solid black;border-bottom: 1px solid gray;border-top: 2px solid gray;">PropVarTempMutPred</th>
<th <th align="center" style="font-weight: normal;border-left: 0px solid black;border-bottom: 1px solid gray;border-top: 2px solid gray;">PropVarTempPhenPred</th>
<th <th align="center" style="font-weight: normal;border-left: 0px solid black;border-bottom: 1px solid gray;border-top: 2px solid gray;">PropVarEnv2MutPred</th>
<th <th align="center" style="font-weight: normal;border-left: 0px solid black;border-right:0px solid black;border-bottom: 1px solid gray;border-top: 2px solid gray;">PropVarEnv2PhenPred</th>
</tr>
<tr>
<td  style="border-left: 0px solid black; ">subsamp_corr_phen_temp          </td>
<td align="right" style="border-left: 0px solid black;background-color: #FFF5F0;">0.02</td>
<td align="right" style="border-left: 0px solid black;background-color: #67000D;color: #FFFFFF;">0.93</td>
<td align="right" style="border-left: 0px solid black;background-color: #FFF5F0;">0.01</td>
<td align="right" style="border-left: 0px solid black;border-right:0px solid black;background-color: #EF3B2C;color: #FFFFFF;">0.54</td>
</tr>
<tr>
<td  style="border-left: 0px solid black; border-top: hidden;">demog_level                     </td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #FEE0D2;">0.06</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #FFF5F0;">0.00</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #FEE0D2;">0.10</td>
<td align="right" style="border-left: 0px solid black;border-right:0px solid black;border-top: hidden;background-color: #FFF5F0;">0.03</td>
</tr>
<tr>
<td  style="border-left: 0px solid black; border-top: hidden;">demog_level_sub                 </td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #FFF5F0;">0.01</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #FFF5F0;">0.02</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #FFF5F0;">0.05</td>
<td align="right" style="border-left: 0px solid black;border-right:0px solid black;border-top: hidden;background-color: #FCBBA1;">0.22</td>
</tr>
<tr>
<td  style="border-left: 0px solid black; border-top: hidden;">arch                            </td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #EF3B2C;color: #FFFFFF;">0.63</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #FFF5F0;">0.00</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #FC9272;">0.34</td>
<td align="right" style="border-left: 0px solid black;border-right:0px solid black;border-top: hidden;background-color: #FFF5F0;">0.00</td>
</tr>
<tr>
<td  style="border-left: 0px solid black; border-top: hidden;">demog_level:demog_level_sub     </td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #FFF5F0;">0.00</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #FFF5F0;">0.00</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #FFF5F0;">0.02</td>
<td align="right" style="border-left: 0px solid black;border-right:0px solid black;border-top: hidden;background-color: #FEE0D2;">0.08</td>
</tr>
<tr>
<td  style="border-left: 0px solid black; border-top: hidden;">demog_level:arch                </td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #FFF5F0;">0.01</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #FFF5F0;">0.00</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #FEE0D2;">0.09</td>
<td align="right" style="border-left: 0px solid black;border-right:0px solid black;border-top: hidden;background-color: #FFF5F0;">0.00</td>
</tr>
<tr>
<td  style="border-left: 0px solid black; border-top: hidden;">demog_level_sub:arch            </td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #FFF5F0;">0.01</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #FFF5F0;">0.00</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #FFF5F0;">0.02</td>
<td align="right" style="border-left: 0px solid black;border-right:0px solid black;border-top: hidden;background-color: #FFF5F0;">0.01</td>
</tr>
<tr>
<td  style="border-left: 0px solid black; border-top: hidden;">demog_level:demog_level_sub:arch</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #FFF5F0;">0.02</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #FFF5F0;">0.00</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #FFF5F0;">0.02</td>
<td align="right" style="border-left: 0px solid black;border-right:0px solid black;border-top: hidden;background-color: #FFF5F0;">0.01</td>
</tr>
<tr>
<td  style="border-left: 0px solid black; border-top: hidden;">Residuals                       </td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #FCBBA1;">0.24</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #FFF5F0;">0.03</td>
<td align="right" style="border-left: 0px solid black;border-top: hidden;background-color: #FC9272;">0.34</td>
<td align="right" style="border-left: 0px solid black;border-right:0px solid black;border-top: hidden;background-color: #FEE0D2;">0.11</td>
</tr>
<tr>
<td colspan="5" align="left" style="font-size:9pt ;border-top: 1px solid black; border-bottom: hidden;"></td>
</tr>
</table>
