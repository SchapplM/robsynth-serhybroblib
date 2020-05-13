% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% palh1m2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [13x1]
%   Generalized joint coordinates (joint angles)
% pkin [20x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi312,phi413,phi710,phi711]';
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
% Datum: 2020-05-02 23:30
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = palh1m2OL_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [13 1]), ...
  'palh1m2OL_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [13x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [20 1]), ...
  'palh1m2OL_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [20x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 21:16:37
% EndTime: 2020-05-02 21:16:38
% DurationCPUTime: 0.38s
% Computational Cost: add. (422->123), mult. (178->104), div. (0->0), fcn. (307->34), ass. (0->81)
t60 = qJ(2) + qJ(7);
t47 = pkin(19) - t60;
t33 = -qJ(10) + t47;
t24 = sin(t33);
t65 = sin(qJ(1));
t88 = t65 * t24;
t27 = cos(t33);
t87 = t65 * t27;
t61 = qJ(2) + qJ(3);
t56 = qJ(4) + t61;
t43 = sin(t56);
t86 = t65 * t43;
t62 = sin(qJ(6));
t85 = t65 * t62;
t63 = sin(qJ(5));
t84 = t65 * t63;
t67 = cos(qJ(5));
t83 = t65 * t67;
t69 = cos(qJ(1));
t82 = t69 * t24;
t81 = t69 * t27;
t80 = t69 * t43;
t79 = t69 * t62;
t78 = t69 * t63;
t77 = t69 * t67;
t59 = qJ(2) + qJ(8);
t64 = sin(qJ(2));
t38 = -t64 * pkin(1) + pkin(15);
t54 = cos(t61);
t9 = pkin(5) * t54 + t38;
t76 = t65 * t9 + 0;
t75 = t69 * t9 + 0;
t58 = pkin(13) + 0;
t49 = sin(t59);
t19 = -pkin(2) * t49 + pkin(15);
t48 = pkin(17) + t61;
t46 = pkin(20) + t60;
t52 = cos(t59);
t74 = pkin(2) * t52 + t58;
t68 = cos(qJ(2));
t18 = t68 * pkin(1) + t58;
t73 = -pkin(16) + t58;
t7 = pkin(4) * sin(t47) + t38;
t72 = pkin(4) * cos(t47) + t18;
t51 = sin(t61);
t71 = pkin(5) * t51 + t18;
t37 = pkin(18) + t56;
t45 = cos(t56);
t70 = pkin(9) * t45 + pkin(11) * t43;
t66 = cos(qJ(6));
t55 = qJ(9) + t59;
t53 = cos(t60);
t50 = sin(t60);
t44 = cos(t55);
t42 = sin(t55);
t41 = t69 * t66;
t40 = t65 * t66;
t39 = t66 * pkin(7) - pkin(14);
t32 = qJ(12) + t48;
t31 = qJ(11) + t46;
t30 = t69 * pkin(15) + 0;
t29 = t65 * pkin(15) + 0;
t28 = qJ(13) + t37;
t26 = cos(t32);
t25 = cos(t31);
t23 = sin(t32);
t22 = sin(t31);
t17 = cos(t28);
t16 = sin(t28);
t15 = t69 * t44;
t14 = t69 * t42;
t13 = t65 * t44;
t12 = t65 * t42;
t11 = t69 * t38 + 0;
t10 = t65 * t38 + 0;
t8 = -pkin(3) * sin(t46) + t38;
t6 = pkin(6) * cos(t48) + t38;
t5 = pkin(12) * t42 + t19;
t2 = pkin(10) * cos(t37) + t9;
t1 = -pkin(8) * t24 + t7;
t3 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t69, -t65, 0, 0; t65, t69, 0, 0; 0, 0, 1, t58; 0, 0, 0, 1; -t69 * t64, -t69 * t68, t65, t30; -t65 * t64, -t65 * t68, -t69, t29; t68, -t64, 0, t58; 0, 0, 0, 1; t69 * t54, -t69 * t51, t65, t11; t65 * t54, -t65 * t51, -t69, t10; t51, t54, 0, t18; 0, 0, 0, 1; t69 * t45, -t80, t65, t75; t65 * t45, -t86, -t69, t76; t43, t45, 0, t71; 0, 0, 0, 1; t45 * t77 + t84, -t45 * t78 + t83, t80, t70 * t69 + t75; t45 * t83 - t78, -t45 * t84 - t77, t86, t70 * t65 + t76; t43 * t67, -t43 * t63, -t45, t43 * pkin(9) - t45 * pkin(11) + t71; 0, 0, 0, 1; t41, -t79, t65, -t69 * pkin(14) + 0; t40, -t85, -t69, -t65 * pkin(14) + 0; t62, t66, 0, t73; 0, 0, 0, 1; -t69 * t50, -t69 * t53, t65, t11; -t65 * t50, -t65 * t53, -t69, t10; t53, -t50, 0, t18; 0, 0, 0, 1; -t69 * t49, -t69 * t52, t65, t30; -t65 * t49, -t65 * t52, -t69, t29; t52, -t49, 0, t58; 0, 0, 0, 1; t14, t15, t65, t69 * t19 + 0; t12, t13, -t69, t65 * t19 + 0; -t44, t42, 0, t74; 0, 0, 0, 1; -t82, t81, t65, t69 * t7 + 0; -t88, t87, -t69, t65 * t7 + 0; -t27, -t24, 0, t72; 0, 0, 0, 1; t69 * t22, t69 * t25, t65, t69 * t8 + 0; t65 * t22, t65 * t25, -t69, t65 * t8 + 0; -t25, t22, 0, pkin(3) * cos(t46) + t18; 0, 0, 0, 1; -t69 * t26, t69 * t23, t65, t69 * t6 + 0; -t65 * t26, t65 * t23, -t69, t65 * t6 + 0; -t23, -t26, 0, pkin(6) * sin(t48) + t18; 0, 0, 0, 1; t69 * t17, -t69 * t16, t65, t69 * t2 + 0; t65 * t17, -t65 * t16, -t69, t65 * t2 + 0; t16, t17, 0, pkin(10) * sin(t37) + t71; 0, 0, 0, 1; t41, -t79, t65, t69 * t39 + 0; t40, -t85, -t69, t65 * t39 + 0; t62, t66, 0, t62 * pkin(7) + t73; 0, 0, 0, 1; t14, t15, t65, t69 * t5 + 0; t12, t13, -t69, t65 * t5 + 0; -t44, t42, 0, -pkin(12) * t44 + t74; 0, 0, 0, 1; t82, -t81, t65, t69 * t1 + 0; t88, -t87, -t69, t65 * t1 + 0; t27, t24, 0, -pkin(8) * t27 + t72; 0, 0, 0, 1;];
T_ges = t3;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,16+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,16+1]); end % symbolisch
for i = 1:16+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
