% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% fourbar2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2]';
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
% Datum: 2020-04-24 20:32
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = fourbar2OL_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbar2OL_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'fourbar2OL_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 20:32:12
% EndTime: 2020-04-24 20:32:12
% DurationCPUTime: 0.09s
% Computational Cost: add. (33->18), mult. (8->6), div. (0->0), fcn. (28->8), ass. (0->14)
t8 = qJ(1) + qJ(2);
t15 = pkin(1) + 0;
t10 = sin(qJ(1));
t14 = t10 * pkin(2) + 0;
t12 = cos(qJ(1));
t13 = t12 * pkin(2) + 0;
t11 = cos(qJ(3));
t9 = sin(qJ(3));
t5 = qJ(4) + t8;
t4 = cos(t8);
t3 = sin(t8);
t2 = cos(t5);
t1 = sin(t5);
t6 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t12, -t10, 0, 0; t10, t12, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; -t4, t3, 0, t13; -t3, -t4, 0, t14; 0, 0, 1, 0; 0, 0, 0, 1; t11, -t9, 0, t15; t9, t11, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; -t2, t1, 0, -pkin(1) * t4 + t13; -t1, -t2, 0, -pkin(1) * t3 + t14; 0, 0, 1, 0; 0, 0, 0, 1; t11, -t9, 0, t11 * pkin(2) + t15; t9, t11, 0, t9 * pkin(2) + 0; 0, 0, 1, 0; 0, 0, 0, 1;];
T_ges = t6;
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
