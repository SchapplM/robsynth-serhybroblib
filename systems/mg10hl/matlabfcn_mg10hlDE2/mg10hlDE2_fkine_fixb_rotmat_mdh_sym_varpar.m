% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% mg10hlDE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [17x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AC,AE,CG,DC,ED,GK,GP,HP,LW,ML,OT,PM,TA,TE,phi23,phi3,phi34]';
% 
% Output:
% T_c_mdh [4x4x(15+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   13:  mdh base (link 0) -> mdh frame (13-1), link (13-1)
%   ...
%   15+1:  mdh base (link 0) -> mdh frame (15)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-11 13:01
% Revision: 6ae2d958c5b90587a0d08029b131cb7b66342a68 (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = mg10hlDE2_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'mg10hlDE2_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [17 1]), ...
  'mg10hlDE2_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [17x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 1
% StartTime: 2020-04-11 13:00:01
% EndTime: 2020-04-11 13:00:02
% DurationCPUTime: 1.09s
% Computational Cost: add. (24491->254), mult. (25650->202), div. (3339->10), fcn. (12414->42), ass. (0->400)
unknown=NaN(64,4);
t1 = cos(qJ(1));
t2 = sin(qJ(1));
t3 = (pkin(11) + 0);
t4 = sin(qJ(2));
t6 = cos(qJ(2));
t14 = cos(pkin(15));
t16 = sin(pkin(15));
t18 = t14 * t4 + t16 * t6;
t20 = -t18 * pkin(2) + pkin(1);
t21 = pkin(1) * t18;
t23 = 0.2e1 * pkin(2) * t21;
t24 = pkin(1) ^ 2;
t34 = sqrt(-(-t23 + t24 + (pkin(2) + pkin(5) - pkin(4)) * (pkin(2) - pkin(5) + pkin(4))) * (-t23 + t24 + (pkin(2) + pkin(5) + pkin(4)) * (pkin(2) - pkin(5) - pkin(4))));
t38 = t14 * t6 - t16 * t4;
t39 = t38 * pkin(2);
t40 = pkin(2) ^ 2;
t41 = pkin(5) ^ 2;
t42 = pkin(4) ^ 2;
t43 = -t23 + t24 + t40 - t41 + t42;
t46 = 0.1e1 / pkin(4);
t48 = -t23 + t24 + t40;
t49 = 0.1e1 / t48;
t56 = atan2(t49 * t46 * (t34 * t20 + t43 * t39), t49 * t46 * (t43 * t20 - t34 * t39));
t57 = qJ(2) + pkin(15) + t56;
t58 = sin(t57);
t59 = t58 * t1;
t60 = cos(t57);
t61 = t60 * t1;
t62 = qJ(2) + pkin(15);
t63 = sin(t62);
t64 = t63 * pkin(1);
t65 = -t64 + pkin(13);
t68 = t58 * t2;
t69 = t60 * t2;
t72 = cos(t62);
t73 = t72 * pkin(1);
t75 = 1 / pkin(7);
t76 = qJ(6) + pkin(8);
t78 = 0.1e1 / t76 * t75;
t79 = pkin(6) - pkin(7) - pkin(8) - qJ(6);
t80 = pkin(6) - pkin(7) + pkin(8) + qJ(6);
t82 = pkin(6) + pkin(7) - pkin(8) - qJ(6);
t83 = pkin(6) + pkin(7) + pkin(8) + qJ(6);
t86 = sqrt(-t83 * t82 * t80 * t79);
t88 = (pkin(6) ^ 2);
t89 = (pkin(7) ^ 2);
t90 = (pkin(8) ^ 2);
t92 = 2 * pkin(8) * qJ(6);
t93 = qJ(6) ^ 2;
t94 = t88 - t89 - t90 - t92 - t93;
t96 = atan2(t86 * t78, t94 * t78);
t97 = qJ(2) + pkin(15) + t56 + pkin(17) - t96 - pkin(16);
t98 = sin(t97);
t100 = cos(t97);
t102 = pkin(15) + t56 + pkin(17) + qJ(2);
t103 = sin(t102);
t104 = t103 * pkin(3);
t105 = -t64 + t104 + pkin(13);
t107 = t105 * t1 + 0;
t111 = t105 * t2 + 0;
t112 = cos(t102);
t113 = t112 * pkin(3);
t114 = t73 - t113 + pkin(11) + 0;
t116 = t76 ^ 2;
t118 = 0.1e1 / t116 / pkin(6);
t123 = t75 * (t88 - t89 + t90 + t92 + t93);
t126 = -t94 * t75 * t86 * t118 + t86 * t123 * t118;
t127 = cos(pkin(16));
t136 = -t75 * t83 * t82 * t80 * t79 * t118 + t94 * t123 * t118;
t137 = sin(pkin(16));
t143 = atan2(t127 * t126 / 0.4e1 + t137 * t136 / 0.4e1, -t127 * t136 / 0.4e1 + t137 * t126 / 0.4e1);
t144 = qJ(2) + pkin(15) + t56 + pkin(17) - t96 - pkin(16) + t143;
t145 = sin(t144);
t146 = t145 * t1;
t147 = cos(t144);
t148 = t147 * t1;
t150 = -t98 * pkin(7) + pkin(13) + t104 - t64;
t151 = t150 * t1;
t153 = t145 * t2;
t154 = t147 * t2;
t155 = t150 * t2;
t157 = t100 * pkin(7);
t159 = cos(qJ(3));
t161 = sin(qJ(3));
t163 = t159 * t148 + t161 * t2;
t166 = -t161 * t148 + t159 * t2;
t167 = pkin(12) * t146;
t171 = -t161 * t1 + t159 * t154;
t174 = -t159 * t1 - t161 * t154;
t175 = pkin(12) * t153;
t177 = t159 * t145;
t178 = t161 * t145;
t179 = pkin(12) * t147;
t181 = sin(qJ(4));
t183 = cos(qJ(4));
t185 = t183 * t146 + t181 * t163;
t188 = -t181 * t146 + t183 * t163;
t189 = pkin(10) * t146;
t193 = t183 * t153 + t181 * t171;
t196 = -t181 * t153 + t183 * t171;
t197 = pkin(10) * t153;
t201 = -t183 * t147 + t181 * t177;
t204 = t181 * t147 + t183 * t177;
t205 = pkin(10) * t147;
t207 = cos(qJ(5));
t209 = sin(qJ(5));
t233 = qJ(2) + pkin(15) + t56 + pkin(17) - pkin(16);
t234 = sin(t233);
t236 = cos(t233);
t242 = t21 - pkin(2);
t243 = -t23 + t24 + t40 + t41 - t42;
t245 = -pkin(1) * t34 * t38 - t243 * t242;
t247 = 0.1e1 / pkin(5);
t252 = t243 * t38 * pkin(1) - t34 * t242;
t253 = t252 ^ 2;
t254 = 0.1e1 / t41;
t256 = t48 ^ 2;
t257 = 0.1e1 / t256;
t259 = t245 ^ 2;
t263 = sqrt(t257 * t254 * t253 + t257 * t254 * t259);
t264 = 0.1e1 / t263;
t265 = t264 * t49 * t247;
t277 = t247 * t252;
t278 = t264 * t49;
t280 = t247 * t245;
t287 = atan2(t34 * t46 * t247, t46 * t247 * (t23 - t24 - t40 + t41 + t42));
t288 = qJ(2) + pkin(15) + t56 + pkin(17) - pkin(16) + t287;
t289 = sin(t288);
t290 = t289 * t1;
t291 = cos(t288);
t292 = t291 * t1;
t294 = -t234 * pkin(6) + pkin(13) + t104 - t64;
t295 = t294 * t1;
t297 = t289 * t2;
t298 = t291 * t2;
t299 = t294 * t2;
t301 = t236 * pkin(6);
t305 = atan2(t49 * t277, t49 * t280);
t306 = atan2(t126 / 0.4e1, t136 / 0.4e1);
t307 = t305 + t306;
t308 = cos(t307);
t310 = sin(t307);
t314 = t264 * t49 * t245 - pkin(14);
t331 = -t58 * pkin(4) + pkin(13) - t64;
unknown(1,1) = 1;
unknown(1,2) = 0;
unknown(1,3) = 0;
unknown(1,4) = 0;
unknown(2,1) = 0;
unknown(2,2) = 1;
unknown(2,3) = 0;
unknown(2,4) = 0;
unknown(3,1) = 0;
unknown(3,2) = 0;
unknown(3,3) = 1;
unknown(3,4) = 0;
unknown(4,1) = 0;
unknown(4,2) = 0;
unknown(4,3) = 0;
unknown(4,4) = 1;
unknown(5,1) = t1;
unknown(5,2) = -t2;
unknown(5,3) = 0;
unknown(5,4) = 0;
unknown(6,1) = t2;
unknown(6,2) = t1;
unknown(6,3) = 0;
unknown(6,4) = 0;
unknown(7,1) = 0;
unknown(7,2) = 0;
unknown(7,3) = 1;
unknown(7,4) = t3;
unknown(8,1) = 0;
unknown(8,2) = 0;
unknown(8,3) = 0;
unknown(8,4) = 1;
unknown(9,1) = -(t4 * t1);
unknown(9,2) = -(t6 * t1);
unknown(9,3) = t2;
unknown(9,4) = (pkin(13) * t1 + 0);
unknown(10,1) = -(t4 * t2);
unknown(10,2) = -(t6 * t2);
unknown(10,3) = -t1;
unknown(10,4) = (pkin(13) * t2 + 0);
unknown(11,1) = t6;
unknown(11,2) = -t4;
unknown(11,3) = 0;
unknown(11,4) = t3;
unknown(12,1) = 0;
unknown(12,2) = 0;
unknown(12,3) = 0;
unknown(12,4) = 1;
unknown(13,1) = t59;
unknown(13,2) = t61;
unknown(13,3) = t2;
unknown(13,4) = (t65 * t1 + 0);
unknown(14,1) = t68;
unknown(14,2) = t69;
unknown(14,3) = -t1;
unknown(14,4) = (t65 * t2 + 0);
unknown(15,1) = -t60;
unknown(15,2) = t58;
unknown(15,3) = 0;
unknown(15,4) = (t73 + pkin(11) + 0);
unknown(16,1) = 0;
unknown(16,2) = 0;
unknown(16,3) = 0;
unknown(16,4) = 1;
unknown(17,1) = -(t98 * t1);
unknown(17,2) = -(t100 * t1);
unknown(17,3) = t2;
unknown(17,4) = t107;
unknown(18,1) = -(t98 * t2);
unknown(18,2) = -(t100 * t2);
unknown(18,3) = -t1;
unknown(18,4) = t111;
unknown(19,1) = t100;
unknown(19,2) = -t98;
unknown(19,3) = 0;
unknown(19,4) = t114;
unknown(20,1) = 0;
unknown(20,2) = 0;
unknown(20,3) = 0;
unknown(20,4) = 1;
unknown(21,1) = t146;
unknown(21,2) = t148;
unknown(21,3) = t2;
unknown(21,4) = (t151 + 0);
unknown(22,1) = t153;
unknown(22,2) = t154;
unknown(22,3) = -t1;
unknown(22,4) = (t155 + 0);
unknown(23,1) = -t147;
unknown(23,2) = t145;
unknown(23,3) = 0;
unknown(23,4) = (t73 - t113 + t157 + pkin(11) + 0);
unknown(24,1) = 0;
unknown(24,2) = 0;
unknown(24,3) = 0;
unknown(24,4) = 1;
unknown(25,1) = t163;
unknown(25,2) = t166;
unknown(25,3) = t146;
unknown(25,4) = (t167 + t151 + 0);
unknown(26,1) = t171;
unknown(26,2) = t174;
unknown(26,3) = t153;
unknown(26,4) = (t175 + t155 + 0);
unknown(27,1) = t177;
unknown(27,2) = -t178;
unknown(27,3) = -t147;
unknown(27,4) = (-t179 + t73 - t113 + t157 + pkin(11) + 0);
unknown(28,1) = 0;
unknown(28,2) = 0;
unknown(28,3) = 0;
unknown(28,4) = 1;
unknown(29,1) = t185;
unknown(29,2) = t188;
unknown(29,3) = t166;
unknown(29,4) = (t189 + t167 + t151 + 0);
unknown(30,1) = t193;
unknown(30,2) = t196;
unknown(30,3) = t174;
unknown(30,4) = (t197 + t175 + t155 + 0);
unknown(31,1) = t201;
unknown(31,2) = t204;
unknown(31,3) = -t178;
unknown(31,4) = (-t205 - t179 + t73 - t113 + t157 + pkin(11) + 0);
unknown(32,1) = 0;
unknown(32,2) = 0;
unknown(32,3) = 0;
unknown(32,4) = 1;
unknown(33,1) = (t209 * t166 + t207 * t188);
unknown(33,2) = (t207 * t166 - t209 * t188);
unknown(33,3) = t185;
unknown(33,4) = (pkin(9) * t185 + t151 + t167 + t189 + 0);
unknown(34,1) = (t209 * t174 + t207 * t196);
unknown(34,2) = (t207 * t174 - t209 * t196);
unknown(34,3) = t193;
unknown(34,4) = (pkin(9) * t193 + t155 + t175 + t197 + 0);
unknown(35,1) = (-t209 * t178 + t207 * t204);
unknown(35,2) = (-t207 * t178 - t209 * t204);
unknown(35,3) = t201;
unknown(35,4) = (pkin(9) * t201 + pkin(11) - t113 + t157 - t179 - t205 + t73 + 0);
unknown(36,1) = 0;
unknown(36,2) = 0;
unknown(36,3) = 0;
unknown(36,4) = 1;
unknown(37,1) = -(t234 * t1);
unknown(37,2) = -(t236 * t1);
unknown(37,3) = t2;
unknown(37,4) = t107;
unknown(38,1) = -(t234 * t2);
unknown(38,2) = -(t236 * t2);
unknown(38,3) = -t1;
unknown(38,4) = t111;
unknown(39,1) = t236;
unknown(39,2) = -t234;
unknown(39,3) = 0;
unknown(39,4) = t114;
unknown(40,1) = 0;
unknown(40,2) = 0;
unknown(40,3) = 0;
unknown(40,4) = 1;
unknown(41,1) = (t265 * t245 * t1);
unknown(41,2) = -(t265 * t252 * t1);
unknown(41,3) = t2;
unknown(41,4) = (-pkin(14) * t1 + 0);
unknown(42,1) = (t265 * t245 * t2);
unknown(42,2) = -(t265 * t252 * t2);
unknown(42,3) = -t1;
unknown(42,4) = (-pkin(14) * t2 + 0);
unknown(43,1) = (t278 * t277);
unknown(43,2) = (t278 * t280);
unknown(43,3) = 0;
unknown(43,4) = t3;
unknown(44,1) = 0;
unknown(44,2) = 0;
unknown(44,3) = 0;
unknown(44,4) = 1;
unknown(45,1) = t290;
unknown(45,2) = t292;
unknown(45,3) = t2;
unknown(45,4) = (t295 + 0);
unknown(46,1) = t297;
unknown(46,2) = t298;
unknown(46,3) = -t1;
unknown(46,4) = (t299 + 0);
unknown(47,1) = -t291;
unknown(47,2) = t289;
unknown(47,3) = 0;
unknown(47,4) = (t73 - t113 + t301 + pkin(11) + 0);
unknown(48,1) = 0;
unknown(48,2) = 0;
unknown(48,3) = 0;
unknown(48,4) = 1;
unknown(49,1) = -(t308 * t1);
unknown(49,2) = (t310 * t1);
unknown(49,3) = t2;
unknown(49,4) = (t314 * t1 + 0);
unknown(50,1) = -(t308 * t2);
unknown(50,2) = (t310 * t2);
unknown(50,3) = -t1;
unknown(50,4) = (t314 * t2 + 0);
unknown(51,1) = -t310;
unknown(51,2) = -t308;
unknown(51,3) = 0;
unknown(51,4) = (t264 * t49 * t252 + pkin(11) + 0);
unknown(52,1) = 0;
unknown(52,2) = 0;
unknown(52,3) = 0;
unknown(52,4) = 1;
unknown(53,1) = -t292;
unknown(53,2) = -t2;
unknown(53,3) = t290;
unknown(53,4) = (qJ(6) * t290 + t295 + 0);
unknown(54,1) = -t298;
unknown(54,2) = t1;
unknown(54,3) = t297;
unknown(54,4) = (qJ(6) * t297 + t299 + 0);
unknown(55,1) = -t289;
unknown(55,2) = 0;
unknown(55,3) = -t291;
unknown(55,4) = (-qJ(6) * t291 + pkin(11) - t113 + t301 + t73 + 0);
unknown(56,1) = 0;
unknown(56,2) = 0;
unknown(56,3) = 0;
unknown(56,4) = 1;
unknown(57,1) = t59;
unknown(57,2) = t61;
unknown(57,3) = t2;
unknown(57,4) = (t331 * t1 + 0);
unknown(58,1) = t68;
unknown(58,2) = t69;
unknown(58,3) = -t1;
unknown(58,4) = (t331 * t2 + 0);
unknown(59,1) = -t60;
unknown(59,2) = t58;
unknown(59,3) = 0;
unknown(59,4) = (t60 * pkin(4) + pkin(11) + t73 + 0);
unknown(60,1) = 0;
unknown(60,2) = 0;
unknown(60,3) = 0;
unknown(60,4) = 1;
unknown(61,1) = -t148;
unknown(61,2) = -t2;
unknown(61,3) = t146;
unknown(61,4) = (-pkin(8) * t146 + t151 + 0);
unknown(62,1) = -t154;
unknown(62,2) = t1;
unknown(62,3) = t153;
unknown(62,4) = (-pkin(8) * t153 + t155 + 0);
unknown(63,1) = -t145;
unknown(63,2) = 0;
unknown(63,3) = -t147;
unknown(63,4) = (pkin(8) * t147 + pkin(11) - t113 + t157 + t73 + 0);
unknown(64,1) = 0;
unknown(64,2) = 0;
unknown(64,3) = 0;
unknown(64,4) = 1;
T_ges = unknown;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,15+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,15+1]); end % symbolisch
for i = 1:15+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
Tc_stack = NaN(3*size(T_c_mdh,3),4);
for i = 1:size(T_c_mdh,3), Tc_stack((i-1)*3+1:3*i,1:4) = T_c_mdh(1:3,1:4,i); end
