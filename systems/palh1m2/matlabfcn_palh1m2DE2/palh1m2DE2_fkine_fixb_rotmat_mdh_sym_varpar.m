% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% palh1m2DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [22x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% 
% Output:
% T_c_mdh [4x4x(16+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   11:  mdh base (link 0) -> mdh frame (11-1), link (11-1)
%   ...
%   16+1:  mdh base (link 0) -> mdh frame (16)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 21:08
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = palh1m2DE2_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2DE2_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2DE2_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [22x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 20:57:31
% EndTime: 2020-05-02 20:57:32
% DurationCPUTime: 1.01s
% Computational Cost: add. (2006->173), mult. (2909->175), div. (0->0), fcn. (4749->70), ass. (0->129)
t113 = sin(qJ(2));
t118 = cos(qJ(2));
t112 = sin(qJ(3));
t117 = cos(qJ(3));
t106 = sin(pkin(20));
t109 = cos(pkin(20));
t115 = sin(pkin(18));
t120 = cos(pkin(18));
t56 = -t120 * t106 + t115 * t109;
t60 = t115 * t106 + t120 * t109;
t122 = t112 * t60 - t56 * t117;
t40 = t112 * t56 + t60 * t117;
t154 = t113 * t40 + t122 * t118;
t121 = cos(pkin(17));
t148 = sin(pkin(17));
t61 = t115 * t121 - t120 * t148;
t64 = t148 * t115 + t121 * t120;
t127 = -t113 * t64 + t61 * t118;
t153 = -t113 * t122 + t40 * t118;
t99 = pkin(22) + pkin(21);
t81 = pkin(18) - pkin(20) - t99;
t74 = qJ(1) + t81;
t152 = -sin(t74) / 0.2e1;
t75 = -qJ(1) + t81;
t151 = -cos(t75) / 0.2e1;
t105 = sin(pkin(22));
t108 = cos(pkin(22));
t58 = t115 * t105 + t108 * t120;
t59 = t120 * t105 - t115 * t108;
t135 = qJ(2) + atan2(-t59 * t113 + t58 * t118, t113 * t58 + t59 * t118);
t24 = pkin(21) - t135;
t150 = pkin(4) * sin(t24);
t102 = qJ(2) + qJ(3);
t93 = cos(t102);
t149 = pkin(5) * t93;
t114 = sin(qJ(1));
t87 = sin(t99);
t88 = cos(t99);
t8 = atan2(t87 * t153 + t154 * t88, -t153 * t88 + t154 * t87) + t102;
t6 = sin(t8);
t147 = t114 * t6;
t119 = cos(qJ(1));
t146 = t119 * t6;
t62 = -t112 * t120 + t117 * t115;
t63 = t112 * t115 + t117 * t120;
t44 = t113 * t63 - t62 * t118;
t45 = t113 * t62 + t63 * t118;
t13 = atan2(t44 * t88 + t87 * t45, -t87 * t44 + t45 * t88);
t12 = -t13 + t24;
t10 = sin(t12);
t142 = t114 * t10;
t11 = cos(t12);
t141 = t114 * t11;
t46 = t113 * t61 + t64 * t118;
t140 = t114 * t46;
t139 = t119 * t10;
t138 = t119 * t11;
t137 = t119 * t46;
t107 = sin(pkin(19));
t110 = cos(pkin(19));
t55 = t117 * t107 + t112 * t110;
t57 = t112 * t107 - t117 * t110;
t50 = qJ(2) + atan2(t57, t55);
t111 = sin(qJ(4));
t134 = t114 * t111;
t116 = cos(qJ(4));
t133 = t114 * t116;
t132 = t119 * t111;
t131 = t119 * t116;
t101 = pkin(18) - pkin(22);
t100 = pkin(13) + 0;
t85 = -t113 * pkin(1) + pkin(15);
t66 = t114 * t85 + 0;
t67 = t119 * t85 + 0;
t79 = t114 * pkin(15) + 0;
t80 = t119 * pkin(15) + 0;
t48 = sin(t50);
t36 = -pkin(2) * t48 + pkin(15);
t23 = pkin(22) + t135;
t89 = pkin(19) + t102;
t130 = t114 * t150 + t66;
t129 = t119 * t150 + t67;
t49 = cos(t50);
t128 = pkin(2) * t49 + t100;
t76 = t118 * pkin(1) + t100;
t126 = -pkin(16) + t100;
t5 = pkin(20) + t8;
t125 = pkin(4) * cos(t24) + t76;
t92 = sin(t102);
t124 = pkin(5) * t92 + t76;
t123 = pkin(10) * cos(t5) + t149;
t104 = qJ(1) - qJ(2);
t103 = qJ(1) + qJ(2);
t95 = qJ(1) - t102;
t94 = qJ(1) + t102;
t91 = -qJ(1) + t101;
t90 = qJ(1) + t101;
t83 = cos(t91);
t82 = sin(t90);
t78 = -cos(t90) / 0.2e1;
t77 = sin(t91) / 0.2e1;
t70 = cos(t74);
t69 = sin(t75);
t65 = t85 + t149;
t54 = pkin(6) * cos(t89) + t85;
t52 = atan2(t57, -t55);
t42 = t52 + t89;
t38 = cos(t42);
t37 = sin(t42);
t35 = t127 * pkin(7) - pkin(14);
t34 = t127 * t114;
t33 = t119 * t127;
t32 = t52 + t50;
t31 = cos(t32);
t30 = sin(t32);
t29 = t119 * t31;
t28 = t119 * t30;
t27 = t114 * t31;
t26 = t114 * t30;
t18 = -pkin(3) * sin(t23) + t85;
t17 = pkin(12) * t30 + t36;
t16 = atan2(t113 * t115 + t118 * t120, -t113 * t120 + t118 * t115) + t23;
t15 = cos(t16);
t14 = sin(t16);
t7 = cos(t8);
t3 = t13 + t5;
t2 = cos(t3);
t1 = sin(t3);
t4 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t119, -t114, 0, 0; t114, t119, 0, 0; 0, 0, 1, t100; 0, 0, 0, 1; -t119 * t113, -t119 * t118, t114, t80; -t114 * t113, -t114 * t118, -t119, t79; t118, -t113, 0, t100; 0, 0, 0, 1; t119 * t93, -t119 * t92, t114, t67; t114 * t93, -t114 * t92, -t119, t66; t92, t93, 0, t76; 0, 0, 0, 1; t119 * t7, -t146, t114, t119 * t65 + 0; t114 * t7, -t147, -t119, t114 * t65 + 0; t6, t7, 0, t124; 0, 0, 0, 1; t7 * t131 + t134, -t7 * t132 + t133, t146, (-t69 / 0.2e1 + t152) * pkin(11) + (-t70 / 0.2e1 + t151) * pkin(9) + (cos(t95) / 0.2e1 + cos(t94) / 0.2e1) * pkin(5) + (-sin(t103) / 0.2e1 + sin(t104) / 0.2e1) * pkin(1) + t80; t7 * t133 - t132, -t7 * t134 - t131, t147, (t70 / 0.2e1 + t151) * pkin(11) + (t69 / 0.2e1 + t152) * pkin(9) + (sin(t94) / 0.2e1 + sin(t95) / 0.2e1) * pkin(5) + (-cos(t104) / 0.2e1 + cos(t103) / 0.2e1) * pkin(1) + t79; t6 * t116, -t6 * t111, -t7, pkin(11) * cos(t81) - pkin(9) * sin(t81) + t124; 0, 0, 0, 1; t33, -t137, t114, -t119 * pkin(14) + 0; t34, -t140, -t119, -t114 * pkin(14) + 0; t46, t127, 0, t126; 0, 0, 0, 1; -t83 / 0.2e1 + t78, t82 / 0.2e1 + t77, t114, t67; -t82 / 0.2e1 + t77, t83 / 0.2e1 + t78, -t119, t66; -sin(t101), -cos(t101), 0, t76; 0, 0, 0, 1; -t119 * t48, -t119 * t49, t114, t80; -t114 * t48, -t114 * t49, -t119, t79; t49, -t48, 0, t100; 0, 0, 0, 1; t28, t29, t114, t119 * t36 + 0; t26, t27, -t119, t114 * t36 + 0; -t31, t30, 0, t128; 0, 0, 0, 1; -t139, t138, t114, t129; -t142, t141, -t119, t130; -t11, -t10, 0, t125; 0, 0, 0, 1; t119 * t14, t119 * t15, t114, t119 * t18 + 0; t114 * t14, t114 * t15, -t119, t114 * t18 + 0; -t15, t14, 0, pkin(3) * cos(t23) + t76; 0, 0, 0, 1; -t119 * t38, t119 * t37, t114, t119 * t54 + 0; -t114 * t38, t114 * t37, -t119, t114 * t54 + 0; -t37, -t38, 0, pkin(6) * sin(t89) + t76; 0, 0, 0, 1; t119 * t2, -t119 * t1, t114, t123 * t119 + t67; t114 * t2, -t114 * t1, -t119, t123 * t114 + t66; t1, t2, 0, pkin(10) * sin(t5) + t124; 0, 0, 0, 1; t33, -t137, t114, t35 * t119 + 0; t34, -t140, -t119, t35 * t114 + 0; t46, t127, 0, t46 * pkin(7) + t126; 0, 0, 0, 1; t28, t29, t114, t119 * t17 + 0; t26, t27, -t119, t114 * t17 + 0; -t31, t30, 0, -pkin(12) * t31 + t128; 0, 0, 0, 1; t139, -t138, t114, -pkin(8) * t139 + t129; t142, -t141, -t119, -pkin(8) * t142 + t130; t11, t10, 0, -pkin(8) * t11 + t125; 0, 0, 0, 1;];
T_ges = t4;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,16+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,16+1]); end % symbolisch
for i = 1:16+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
Tc_stack = NaN(3*size(T_c_mdh,3),4);
for i = 1:size(T_c_mdh,3), Tc_stack((i-1)*3+1:3*i,1:4) = T_c_mdh(1:3,1:4,i); end
