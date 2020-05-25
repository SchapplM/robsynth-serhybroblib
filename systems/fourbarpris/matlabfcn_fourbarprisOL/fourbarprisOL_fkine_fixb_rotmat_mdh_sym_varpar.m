% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% fourbarprisOL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[GK,GP,HP]';
% 
% Output:
% T_c_mdh [4x4x(5+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   4:  mdh base (link 0) -> mdh frame (4-1), link (4-1)
%   ...
%   5+1:  mdh base (link 0) -> mdh frame (5)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 09:52
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = fourbarprisOL_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbarprisOL_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisOL_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:51:59
% EndTime: 2020-05-07 09:51:59
% DurationCPUTime: 0.07s
% Computational Cost: add. (27->23), mult. (14->10), div. (0->0), fcn. (38->6), ass. (0->11)
t10 = pkin(1) + 0;
t9 = cos(qJ(1));
t8 = cos(qJ(3));
t7 = cos(qJ(4));
t6 = sin(qJ(1));
t5 = sin(qJ(3));
t4 = sin(qJ(4));
t3 = -pkin(3) - qJ(2);
t2 = -t5 * t4 + t8 * t7;
t1 = t8 * t4 + t5 * t7;
t11 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; -t9, t6, 0, t10; -t6, -t9, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; -t6, 0, -t9, -t9 * qJ(2) + t10; t9, 0, -t6, -t6 * qJ(2) + 0; 0, -1, 0, 0; 0, 0, 0, 1; -t8, t5, 0, 0; -t5, -t8, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t2, -t1, 0, -t8 * pkin(2) + 0; t1, t2, 0, -t5 * pkin(2) + 0; 0, 0, 1, 0; 0, 0, 0, 1; -t9, t6, 0, t3 * t9 + t10; -t6, -t9, 0, t3 * t6 + 0; 0, 0, 1, 0; 0, 0, 0, 1;];
T_ges = t11;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
Tc_stack = NaN(3*size(T_c_mdh,3),4);
for i = 1:size(T_c_mdh,3), Tc_stack((i-1)*3+1:3*i,1:4) = T_c_mdh(1:3,1:4,i); end
