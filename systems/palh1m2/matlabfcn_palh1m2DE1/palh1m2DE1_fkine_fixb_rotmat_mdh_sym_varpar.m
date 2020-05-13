% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% palh1m2DE1
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
% Datum: 2020-05-01 21:04
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = palh1m2DE1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2DE1_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2DE1_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [22x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-01 20:55:33
% EndTime: 2020-05-01 20:55:34
% DurationCPUTime: 0.62s
% Computational Cost: add. (446->165), mult. (637->157), div. (0->0), fcn. (902->24), ass. (0->103)
t119 = sin(pkin(17));
t77 = sin(pkin(18));
t82 = cos(pkin(18));
t83 = cos(pkin(17));
t33 = -t82 * t119 + t77 * t83;
t35 = t119 * t77 + t83 * t82;
t75 = sin(qJ(2));
t80 = cos(qJ(2));
t14 = t75 * t33 + t35 * t80;
t81 = cos(qJ(1));
t127 = t14 * t81;
t76 = sin(qJ(1));
t126 = t76 * t14;
t65 = sin(pkin(22));
t69 = cos(pkin(22));
t125 = t77 * t65 + t69 * t82;
t92 = t33 * t80 - t75 * t35;
t124 = t76 * t92;
t123 = t92 * t81;
t67 = sin(pkin(20));
t71 = cos(pkin(20));
t27 = -t82 * t67 + t77 * t71;
t31 = t77 * t67 + t82 * t71;
t62 = pkin(22) + pkin(21);
t57 = sin(t62);
t58 = cos(t62);
t5 = -t27 * t58 + t57 * t31;
t121 = t5 * t76;
t120 = t5 * t81;
t68 = sin(pkin(19));
t72 = cos(pkin(19));
t74 = sin(qJ(3));
t79 = cos(qJ(3));
t24 = t79 * t68 + t74 * t72;
t118 = t24 * t80;
t28 = -t74 * t68 + t79 * t72;
t117 = t28 * t80;
t115 = t75 * t28;
t113 = t75 * t74;
t104 = t79 * t80;
t32 = t104 - t113;
t112 = t76 * t32;
t105 = t79 * t75;
t34 = t80 * t74 + t105;
t111 = t76 * t34;
t64 = qJ(3) + qJ(2);
t59 = sin(t64);
t110 = t76 * t59;
t60 = cos(t64);
t109 = t76 * t60;
t108 = t76 * t75;
t107 = t76 * t80;
t103 = t81 * t32;
t102 = t81 * t34;
t101 = t81 * t59;
t100 = t81 * t60;
t99 = t81 * t75;
t98 = t81 * t80;
t30 = t82 * t65 - t77 * t69;
t54 = -pkin(1) * t75 + pkin(15);
t66 = sin(pkin(21));
t70 = cos(pkin(21));
t3 = (-t125 * t70 + t30 * t66) * pkin(4) + t54;
t97 = t3 * t76 + 0;
t96 = t3 * t81 + 0;
t63 = pkin(13) + 0;
t56 = t74 * pkin(5) + pkin(1);
t20 = pkin(5) * t104 - t56 * t75 + pkin(15);
t95 = t20 * t76 + 0;
t94 = t20 * t81 + 0;
t93 = t68 * t113;
t91 = pkin(2) * t72 * t105 + t63;
t45 = t80 * pkin(1) + t63;
t90 = -pkin(16) + t63;
t89 = pkin(5) * t105 + t56 * t80 + t63;
t42 = -t77 * pkin(9) + t82 * pkin(11);
t43 = t82 * pkin(9) + t77 * pkin(11);
t16 = t42 * t71 + t67 * t43;
t17 = -t67 * t42 + t43 * t71;
t88 = t16 * t57 - t17 * t58;
t25 = t68 * t75 - t72 * t80;
t29 = t68 * t80 + t72 * t75;
t10 = t25 * t79 + t74 * t29;
t4 = t27 * t57 + t58 * t31;
t87 = pkin(10) * (-t57 * t77 - t58 * t82);
t73 = sin(qJ(4));
t86 = t73 * t4;
t78 = cos(qJ(4));
t85 = t4 * t78;
t84 = t45 + (t125 * t66 + t30 * t70) * pkin(4);
t47 = t81 * pkin(15) + 0;
t46 = t76 * pkin(15) + 0;
t38 = -pkin(3) * t82 + t54;
t37 = t54 * t81 + 0;
t36 = t54 * t76 + 0;
t22 = -t24 * pkin(2) - pkin(12);
t21 = -t24 * pkin(6) - pkin(1);
t11 = -t25 * t74 + t29 * t79;
t9 = t92 * pkin(7) - pkin(14);
t8 = pkin(2) * t117 + t22 * t75 + pkin(15);
t7 = pkin(6) * t117 + t21 * t75 + pkin(15);
t6 = pkin(15) + (-t24 * t75 + t117) * pkin(2);
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t81, -t76, 0, 0; t76, t81, 0, 0; 0, 0, 1, t63; 0, 0, 0, 1; -t99, -t98, t76, t47; -t108, -t107, -t81, t46; t80, -t75, 0, t63; 0, 0, 0, 1; t103, -t102, t76, t37; t112, -t111, -t81, t36; t34, t32, 0, t45; 0, 0, 0, 1; -t4 * t81, -t120, t76, t94; -t4 * t76, -t121, -t81, t95; t5, -t4, 0, t89; 0, 0, 0, 1; t76 * t73 - t81 * t85, t76 * t78 + t81 * t86, t120, t88 * t81 + t94; -t81 * t73 - t76 * t85, t76 * t86 - t81 * t78, t121, t88 * t76 + t95; t5 * t78, -t5 * t73, t4, t16 * t58 + t17 * t57 + t89; 0, 0, 0, 1; t123, -t127, t76, -t81 * pkin(14) + 0; t124, -t126, -t81, -t76 * pkin(14) + 0; t14, t92, 0, t90; 0, 0, 0, 1; -t81 * t125, -t81 * t30, t76, t37; -t76 * t125, -t76 * t30, -t81, t36; t30, -t125, 0, t45; 0, 0, 0, 1; -t81 * t10, -t81 * t11, t76, t47; -t10 * t76, (-t115 - t118) * t76, -t81, t46; t11, -t10, 0, t63; 0, 0, 0, 1; -t99, -t98, t76, t6 * t81 + 0; -t108, -t107, -t81, t6 * t76 + 0; t80, -t75, 0, (-t93 + t118) * pkin(2) + t91; 0, 0, 0, 1; t100, -t101, t76, t96; t109, -t110, -t81, t97; t59, t60, 0, t84; 0, 0, 0, 1; -t99, -t98, t76, t38 * t81 + 0; -t108, -t107, -t81, t38 * t76 + 0; t80, -t75, 0, -pkin(3) * t77 + t45; 0, 0, 0, 1; -t99, -t98, t76, t7 * t81 + 0; -t108, -t107, -t81, t7 * t76 + 0; t80, -t75, 0, pkin(6) * t115 - t21 * t80 + t63; 0, 0, 0, 1; -t103, t102, t76, t81 * t87 + t94; -t112, t111, -t81, t76 * t87 + t95; -t34, -t32, 0, (t57 * t82 - t58 * t77) * pkin(10) + t89; 0, 0, 0, 1; t123, -t127, t76, t9 * t81 + 0; t124, -t126, -t81, t9 * t76 + 0; t14, t92, 0, t14 * pkin(7) + t90; 0, 0, 0, 1; -t99, -t98, t76, t8 * t81 + 0; -t108, -t107, -t81, t8 * t76 + 0; t80, -t75, 0, -pkin(2) * t93 - t22 * t80 + t91; 0, 0, 0, 1; -t100, t101, t76, pkin(8) * t100 + t96; -t109, t110, -t81, pkin(8) * t109 + t97; -t59, -t60, 0, t59 * pkin(8) + t84; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,16+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,16+1]); end % symbolisch
for i = 1:16+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
