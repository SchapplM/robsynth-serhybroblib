% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% fourbar2DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
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
% Datum: 2020-04-24 20:27
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = fourbar2DE2_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar2DE2_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'fourbar2DE2_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 20:27:45
% EndTime: 2020-04-24 20:27:45
% DurationCPUTime: 0.08s
% Computational Cost: add. (27->13), mult. (6->2), div. (0->0), fcn. (46->6), ass. (0->9)
t8 = cos(qJ(1));
t9 = pkin(2) * t8 + 0;
t7 = sin(qJ(1));
t5 = pkin(2) * t7 + 0;
t4 = pkin(1) + t9;
t3 = qJ(1) + atan2(t7, -t8) + atan2(t7, t8);
t2 = cos(t3);
t1 = sin(t3);
t6 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t8, -t7, 0, 0; t7, t8, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; 1, 0, 0, t9; 0, 1, 0, t5; 0, 0, 1, 0; 0, 0, 0, 1; t8, -t7, 0, pkin(1) + 0; t7, t8, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; -t2, t1, 0, t4; -t1, -t2, 0, t5; 0, 0, 1, 0; 0, 0, 0, 1; t8, -t7, 0, t4; t7, t8, 0, t5; 0, 0, 1, 0; 0, 0, 0, 1;];
T_ges = t6;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
