% Calculate homogenous joint transformation matrices for
% palh3m2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [10x1]
%   Generalized joint coordinates (joint angles)
% pkin [16x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi410,phi78,phi79]';
% 
% Output:
% T_mdh [4x4x12]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)
% T_stack [(12+1)*3 x 4]
%   stacked matrices from T_mdh into one 2D array, last row left out.
%   Last row only contains [0 0 0 1].

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-30 21:24
% Revision: b9e8aa5c608190a7b43c48aaebfd2074f0379b0d (2020-06-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_mdh, T_stack] = palh3m2OL_joint_trafo_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m2OL_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [10x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m2OL_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [16x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-30 21:18:45
% EndTime: 2020-06-30 21:18:45
% DurationCPUTime: 0.11s
% Computational Cost: add. (24->21), mult. (30->18), div. (0->0), fcn. (82->26), ass. (0->33)
t98 = cos(qJ(1));
t97 = cos(qJ(2));
t96 = cos(qJ(3));
t95 = cos(qJ(4));
t94 = cos(qJ(5));
t93 = cos(qJ(6));
t92 = cos(qJ(7));
t91 = cos(qJ(8));
t90 = cos(qJ(9));
t89 = sin(qJ(1));
t88 = sin(qJ(2));
t87 = sin(qJ(3));
t86 = sin(qJ(4));
t85 = sin(qJ(5));
t84 = sin(qJ(6));
t83 = sin(qJ(7));
t82 = sin(qJ(8));
t81 = sin(qJ(9));
t80 = cos(pkin(15));
t79 = cos(pkin(16));
t78 = cos(qJ(10));
t77 = sin(pkin(15));
t76 = sin(pkin(16));
t75 = sin(qJ(10));
t74 = cos(pkin(14));
t73 = sin(pkin(14));
t72 = -t77 * t82 - t80 * t91;
t71 = t76 * t81 - t79 * t90;
t70 = t77 * t91 - t80 * t82;
t69 = -t76 * t90 - t79 * t81;
t68 = -t73 * t75 + t74 * t78;
t67 = t73 * t78 + t74 * t75;
t1 = [t98, -t89, 0, 0; t89, t98, 0, 0; 0, 0, 1, pkin(11); t97, -t88, 0, pkin(12); 0, 0, -1, 0; t88, t97, 0, 0; -t96, t87, 0, pkin(1); -t87, -t96, 0, 0; 0, 0, 1, 0; t95, -t86, 0, pkin(4); t86, t95, 0, 0; 0, 0, 1, 0; t94, -t85, 0, pkin(8); 0, 0, -1, -pkin(10); t85, t94, 0, 0; t93, -t84, 0, -pkin(6); 0, 0, -1, 0; t84, t93, 0, pkin(13); t92, -t83, 0, pkin(1); t83, t92, 0, 0; 0, 0, 1, 0; t72, -t70, 0, t80 * pkin(3); t70, t72, 0, -t77 * pkin(3); 0, 0, 1, 0; t71, -t69, 0, t79 * pkin(2); t69, t71, 0, t76 * pkin(2); 0, 0, 1, 0; t68, -t67, 0, t74 * pkin(9); t67, t68, 0, t73 * pkin(9); 0, 0, 1, 0; 1, 0, 0, pkin(5); 0, 1, 0, 0; 0, 0, 1, 0; -1, 0, 0, pkin(7); 0, -1, 0, 0; 0, 0, 1, 0;];
T_stack = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,12);             % numerisch
else,                         T_mdh = sym('xx', [4,4,12]); end % symbolisch

for i = 1:12
  T_mdh(:,:,i) = [T_stack((i-1)*3+1 : 3*i, :);[0 0 0 1]];
end
