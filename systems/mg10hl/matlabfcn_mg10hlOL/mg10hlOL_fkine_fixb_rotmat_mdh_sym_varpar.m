% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% mg10hlOL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [13x1]
%   Generalized joint coordinates (joint angles)
% pkin [16x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AC,CG,DC,ED,GK,GP,HP,LW,ML,OT,PM,TA,TE,phi23,phi3,phi34]';
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
% Datum: 2020-04-11 13:06
% Revision: 6ae2d958c5b90587a0d08029b131cb7b66342a68 (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = mg10hlOL_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'mg10hlOL_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [13x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'mg10hlOL_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [16x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 1
% StartTime: 2020-04-11 13:05:44
% EndTime: 2020-04-11 13:05:44
% DurationCPUTime: 0.36s
% Computational Cost: add. (4184->281), mult. (6926->175), div. (0->0), fcn. (9902->28), ass. (0->373)
unknown=NaN(64,4);
t1 = cos(qJ(1));
t2 = sin(qJ(1));
t3 = (pkin(10) + 0);
t4 = sin(qJ(2));
t5 = t1 * t4;
t6 = cos(qJ(2));
t7 = t1 * t6;
t8 = t1 * pkin(12);
t10 = t2 * t4;
t11 = t2 * t6;
t12 = t2 * pkin(12);
t14 = cos(pkin(14));
t15 = cos(qJ(3));
t17 = sin(pkin(14));
t18 = sin(qJ(3));
t20 = -t14 * t15 + t17 * t18;
t24 = -t14 * t18 - t15 * t17;
t26 = -t20 * t5 - t24 * t7;
t29 = -t20 * t7 + t24 * t5;
t30 = t14 * pkin(1);
t31 = t5 * t30;
t32 = t17 * pkin(1);
t33 = t7 * t32;
t37 = -t10 * t20 - t11 * t24;
t40 = t10 * t24 - t11 * t20;
t41 = t10 * t30;
t42 = t11 * t32;
t46 = t20 * t6 - t24 * t4;
t49 = -t20 * t4 - t24 * t6;
t51 = t6 * t14 * pkin(1);
t53 = t4 * t17 * pkin(1);
t55 = cos(pkin(16));
t56 = -qJ(4) + pkin(15);
t57 = cos(t56);
t59 = sin(pkin(16));
t60 = sin(t56);
t62 = -t55 * t57 - t59 * t60;
t66 = t55 * t60 - t57 * t59;
t68 = t26 * t62 + t29 * t66;
t71 = -t26 * t66 + t29 * t62;
t73 = t26 * t55 * pkin(2);
t75 = t29 * t59 * pkin(2);
t76 = t73 + t75 - t31 - t33 + t8 + 0;
t79 = t37 * t62 + t40 * t66;
t82 = -t37 * t66 + t40 * t62;
t84 = t37 * t55 * pkin(2);
t86 = t40 * t59 * pkin(2);
t87 = t84 + t86 - t41 - t42 + t12 + 0;
t90 = t46 * t62 + t49 * t66;
t93 = -t46 * t66 + t49 * t62;
t95 = t46 * t55 * pkin(2);
t97 = t49 * t59 * pkin(2);
t98 = t95 + t97 + t51 - t53 + pkin(10) + 0;
t99 = cos(qJ(5));
t101 = sin(qJ(5));
t103 = -t101 * t71 - t68 * t99;
t106 = t101 * t68 - t71 * t99;
t107 = t68 * pkin(6);
t111 = -t101 * t82 - t79 * t99;
t114 = t101 * t79 - t82 * t99;
t115 = t79 * pkin(6);
t119 = -t101 * t93 - t90 * t99;
t122 = t101 * t90 - t93 * t99;
t123 = t90 * pkin(6);
t125 = cos(qJ(6));
t127 = sin(qJ(6));
t129 = t106 * t125 + t127 * t2;
t132 = -t106 * t127 + t125 * t2;
t133 = t103 * pkin(11);
t137 = -t1 * t127 + t114 * t125;
t140 = -t1 * t125 - t114 * t127;
t141 = t111 * pkin(11);
t143 = t122 * t125;
t144 = t122 * t127;
t145 = t119 * pkin(11);
t147 = sin(qJ(7));
t149 = cos(qJ(7));
t151 = t103 * t149 + t129 * t147;
t154 = -t103 * t147 + t129 * t149;
t155 = t103 * pkin(9);
t159 = t111 * t149 + t137 * t147;
t162 = -t111 * t147 + t137 * t149;
t163 = t111 * pkin(9);
t167 = t119 * t149 + t143 * t147;
t170 = -t119 * t147 + t143 * t149;
t171 = t119 * pkin(9);
t173 = cos(qJ(8));
t175 = sin(qJ(8));
t199 = cos(qJ(9));
t201 = sin(qJ(9));
t203 = t199 * t68 + t201 * t71;
t206 = t199 * t71 - t201 * t68;
t209 = t199 * t79 + t201 * t82;
t212 = t199 * t82 - t201 * t79;
t215 = t199 * t90 + t201 * t93;
t218 = t199 * t93 - t201 * t90;
t219 = cos(qJ(10));
t220 = t1 * t219;
t221 = sin(qJ(10));
t222 = t1 * t221;
t223 = t1 * pkin(13);
t225 = t2 * t219;
t226 = t2 * t221;
t227 = t2 * pkin(13);
t229 = cos(qJ(11));
t231 = sin(qJ(11));
t233 = -t203 * t229 - t206 * t231;
t236 = t203 * t231 - t206 * t229;
t237 = t203 * pkin(5);
t241 = -t209 * t229 - t212 * t231;
t244 = t209 * t231 - t212 * t229;
t245 = t209 * pkin(5);
t249 = -t215 * t229 - t218 * t231;
t252 = t215 * t231 - t218 * t229;
t253 = t215 * pkin(5);
t255 = cos(qJ(12));
t257 = sin(qJ(12));
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
unknown(9,1) = -t5;
unknown(9,2) = -t7;
unknown(9,3) = t2;
unknown(9,4) = (t8 + 0);
unknown(10,1) = -t10;
unknown(10,2) = -t11;
unknown(10,3) = -t1;
unknown(10,4) = (t12 + 0);
unknown(11,1) = t6;
unknown(11,2) = -t4;
unknown(11,3) = 0;
unknown(11,4) = t3;
unknown(12,1) = 0;
unknown(12,2) = 0;
unknown(12,3) = 0;
unknown(12,4) = 1;
unknown(13,1) = t26;
unknown(13,2) = t29;
unknown(13,3) = t2;
unknown(13,4) = (-t31 - t33 + t8 + 0);
unknown(14,1) = t37;
unknown(14,2) = t40;
unknown(14,3) = -t1;
unknown(14,4) = (-t41 - t42 + t12 + 0);
unknown(15,1) = t46;
unknown(15,2) = t49;
unknown(15,3) = 0;
unknown(15,4) = (t51 - t53 + pkin(10) + 0);
unknown(16,1) = 0;
unknown(16,2) = 0;
unknown(16,3) = 0;
unknown(16,4) = 1;
unknown(17,1) = t68;
unknown(17,2) = t71;
unknown(17,3) = t2;
unknown(17,4) = t76;
unknown(18,1) = t79;
unknown(18,2) = t82;
unknown(18,3) = -t1;
unknown(18,4) = t87;
unknown(19,1) = t90;
unknown(19,2) = t93;
unknown(19,3) = 0;
unknown(19,4) = t98;
unknown(20,1) = 0;
unknown(20,2) = 0;
unknown(20,3) = 0;
unknown(20,4) = 1;
unknown(21,1) = t103;
unknown(21,2) = t106;
unknown(21,3) = t2;
unknown(21,4) = (t107 + t73 + t75 - t31 - t33 + t8 + 0);
unknown(22,1) = t111;
unknown(22,2) = t114;
unknown(22,3) = -t1;
unknown(22,4) = (t115 + t84 + t86 - t41 - t42 + t12 + 0);
unknown(23,1) = t119;
unknown(23,2) = t122;
unknown(23,3) = 0;
unknown(23,4) = (t123 + t95 + t97 + t51 - t53 + pkin(10) + 0);
unknown(24,1) = 0;
unknown(24,2) = 0;
unknown(24,3) = 0;
unknown(24,4) = 1;
unknown(25,1) = t129;
unknown(25,2) = t132;
unknown(25,3) = t103;
unknown(25,4) = (t133 + t107 + t73 + t75 - t31 - t33 + t8 + 0);
unknown(26,1) = t137;
unknown(26,2) = t140;
unknown(26,3) = t111;
unknown(26,4) = (t141 + t115 + t84 + t86 - t41 - t42 + t12 + 0);
unknown(27,1) = t143;
unknown(27,2) = -t144;
unknown(27,3) = t119;
unknown(27,4) = (t145 + t123 + t95 + t97 + t51 - t53 + pkin(10) + 0);
unknown(28,1) = 0;
unknown(28,2) = 0;
unknown(28,3) = 0;
unknown(28,4) = 1;
unknown(29,1) = t151;
unknown(29,2) = t154;
unknown(29,3) = t132;
unknown(29,4) = (t155 + t133 + t107 + t73 + t75 - t31 - t33 + t8 + 0);
unknown(30,1) = t159;
unknown(30,2) = t162;
unknown(30,3) = t140;
unknown(30,4) = (t163 + t141 + t115 + t84 + t86 - t41 - t42 + t12 + 0);
unknown(31,1) = t167;
unknown(31,2) = t170;
unknown(31,3) = -t144;
unknown(31,4) = (t171 + t145 + t123 + t95 + t97 + t51 - t53 + pkin(10) + 0);
unknown(32,1) = 0;
unknown(32,2) = 0;
unknown(32,3) = 0;
unknown(32,4) = 1;
unknown(33,1) = (t132 * t175 + t154 * t173);
unknown(33,2) = (t132 * t173 - t154 * t175);
unknown(33,3) = t151;
unknown(33,4) = (pkin(8) * t151 + t107 + t133 + t155 - t31 - t33 + t73 + t75 + t8 + 0);
unknown(34,1) = (t140 * t175 + t162 * t173);
unknown(34,2) = (t140 * t173 - t162 * t175);
unknown(34,3) = t159;
unknown(34,4) = (pkin(8) * t159 + t115 + t12 + t141 + t163 - t41 - t42 + t84 + t86 + 0);
unknown(35,1) = (-t144 * t175 + t170 * t173);
unknown(35,2) = (-t144 * t173 - t170 * t175);
unknown(35,3) = t167;
unknown(35,4) = (pkin(8) * t167 + pkin(10) + t123 + t145 + t171 + t51 - t53 + t95 + t97 + 0);
unknown(36,1) = 0;
unknown(36,2) = 0;
unknown(36,3) = 0;
unknown(36,4) = 1;
unknown(37,1) = t203;
unknown(37,2) = t206;
unknown(37,3) = t2;
unknown(37,4) = t76;
unknown(38,1) = t209;
unknown(38,2) = t212;
unknown(38,3) = -t1;
unknown(38,4) = t87;
unknown(39,1) = t215;
unknown(39,2) = t218;
unknown(39,3) = 0;
unknown(39,4) = t98;
unknown(40,1) = 0;
unknown(40,2) = 0;
unknown(40,3) = 0;
unknown(40,4) = 1;
unknown(41,1) = t220;
unknown(41,2) = -t222;
unknown(41,3) = t2;
unknown(41,4) = (-t223 + 0);
unknown(42,1) = t225;
unknown(42,2) = -t226;
unknown(42,3) = -t1;
unknown(42,4) = (-t227 + 0);
unknown(43,1) = t221;
unknown(43,2) = t219;
unknown(43,3) = 0;
unknown(43,4) = t3;
unknown(44,1) = 0;
unknown(44,2) = 0;
unknown(44,3) = 0;
unknown(44,4) = 1;
unknown(45,1) = t233;
unknown(45,2) = t236;
unknown(45,3) = t2;
unknown(45,4) = (t237 + t73 + t75 - t31 - t33 + t8 + 0);
unknown(46,1) = t241;
unknown(46,2) = t244;
unknown(46,3) = -t1;
unknown(46,4) = (t245 + t84 + t86 - t41 - t42 + t12 + 0);
unknown(47,1) = t249;
unknown(47,2) = t252;
unknown(47,3) = 0;
unknown(47,4) = (t253 + t95 + t97 + t51 - t53 + pkin(10) + 0);
unknown(48,1) = 0;
unknown(48,2) = 0;
unknown(48,3) = 0;
unknown(48,4) = 1;
unknown(49,1) = (-t220 * t255 + t222 * t257);
unknown(49,2) = (t220 * t257 + t222 * t255);
unknown(49,3) = t2;
unknown(49,4) = (pkin(4) * t220 - t223 + 0);
unknown(50,1) = (-t225 * t255 + t226 * t257);
unknown(50,2) = (t225 * t257 + t226 * t255);
unknown(50,3) = -t1;
unknown(50,4) = (pkin(4) * t225 - t227 + 0);
unknown(51,1) = (-t219 * t257 - t221 * t255);
unknown(51,2) = (-t219 * t255 + t221 * t257);
unknown(51,3) = 0;
unknown(51,4) = (pkin(4) * t221 + pkin(10) + 0);
unknown(52,1) = 0;
unknown(52,2) = 0;
unknown(52,3) = 0;
unknown(52,4) = 1;
unknown(53,1) = -t236;
unknown(53,2) = -t2;
unknown(53,3) = t233;
unknown(53,4) = (qJ(13) * t233 + t237 - t31 - t33 + t73 + t75 + t8 + 0);
unknown(54,1) = -t244;
unknown(54,2) = t1;
unknown(54,3) = t241;
unknown(54,4) = (qJ(13) * t241 + t12 + t245 - t41 - t42 + t84 + t86 + 0);
unknown(55,1) = -t252;
unknown(55,2) = 0;
unknown(55,3) = t249;
unknown(55,4) = (qJ(13) * t249 + pkin(10) + t253 + t51 - t53 + t95 + t97 + 0);
unknown(56,1) = 0;
unknown(56,2) = 0;
unknown(56,3) = 0;
unknown(56,4) = 1;
unknown(57,1) = t26;
unknown(57,2) = t29;
unknown(57,3) = t2;
unknown(57,4) = (-pkin(3) * t26 - t31 - t33 + t8 + 0);
unknown(58,1) = t37;
unknown(58,2) = t40;
unknown(58,3) = -t1;
unknown(58,4) = (-pkin(3) * t37 + t12 - t41 - t42 + 0);
unknown(59,1) = t46;
unknown(59,2) = t49;
unknown(59,3) = 0;
unknown(59,4) = (-pkin(3) * t46 + pkin(10) + t51 - t53 + 0);
unknown(60,1) = 0;
unknown(60,2) = 0;
unknown(60,3) = 0;
unknown(60,4) = 1;
unknown(61,1) = -t106;
unknown(61,2) = -t2;
unknown(61,3) = t103;
unknown(61,4) = (-pkin(7) * t103 + t107 - t31 - t33 + t73 + t75 + t8 + 0);
unknown(62,1) = -t114;
unknown(62,2) = t1;
unknown(62,3) = t111;
unknown(62,4) = (-pkin(7) * t111 + t115 + t12 - t41 - t42 + t84 + t86 + 0);
unknown(63,1) = -t122;
unknown(63,2) = 0;
unknown(63,3) = t119;
unknown(63,4) = (-pkin(7) * t119 + pkin(10) + t123 + t51 - t53 + t95 + t97 + 0);
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
