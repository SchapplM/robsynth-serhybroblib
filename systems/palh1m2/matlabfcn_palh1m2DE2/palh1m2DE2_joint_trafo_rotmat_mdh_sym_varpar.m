% Calculate homogenous joint transformation matrices for
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
% T_mdh [4x4x16]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 21:08
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_mdh = palh1m2DE2_joint_trafo_rotmat_mdh_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2DE2_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2DE2_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [22x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 20:57:32
% EndTime: 2020-05-02 20:57:32
% DurationCPUTime: 0.19s
% Computational Cost: add. (206->43), mult. (336->44), div. (0->0), fcn. (552->24), ass. (0->43)
t183 = sin(pkin(20));
t186 = cos(pkin(20));
t192 = sin(pkin(18));
t198 = cos(pkin(18));
t168 = -t198 * t183 + t192 * t186;
t172 = t192 * t183 + t198 * t186;
t189 = sin(qJ(3));
t195 = cos(qJ(3));
t162 = t168 * t195 - t189 * t172;
t165 = t189 * t168 + t172 * t195;
t190 = sin(qJ(2));
t196 = cos(qJ(2));
t159 = t162 * t196 - t190 * t165;
t160 = t190 * t162 + t165 * t196;
t181 = pkin(22) + pkin(21);
t179 = sin(t181);
t180 = cos(t181);
t155 = t159 * t180 - t179 * t160;
t157 = t159 * t179 + t160 * t180;
t182 = sin(pkin(22));
t185 = cos(pkin(22));
t170 = t192 * t182 + t185 * t198;
t171 = t198 * t182 - t192 * t185;
t164 = t170 * t196 - t190 * t171;
t193 = sin(pkin(17));
t199 = cos(pkin(17));
t173 = t192 * t199 - t198 * t193;
t174 = t193 * t192 + t199 * t198;
t200 = t190 * t173 + t174 * t196;
t197 = cos(qJ(1));
t194 = cos(qJ(4));
t191 = sin(qJ(1));
t188 = sin(qJ(4));
t187 = cos(pkin(19));
t184 = sin(pkin(19));
t178 = pkin(18) - pkin(22) - qJ(2) - qJ(3);
t176 = cos(t178);
t175 = sin(t178);
t169 = -t189 * t184 + t195 * t187;
t167 = t195 * t184 + t189 * t187;
t166 = t173 * t196 - t174 * t190;
t163 = t190 * t170 + t171 * t196;
t1 = [t197, -t191, 0, 0; t191, t197, 0, 0; 0, 0, 1, pkin(13); 0, 0, 0, 1; -t190, -t196, 0, pkin(15); 0, 0, -1, 0; t196, -t190, 0, 0; 0, 0, 0, 1; t189, t195, 0, pkin(1); -t195, t189, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; -t157, t155, 0, pkin(5); -t155, -t157, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t194, -t188, 0, pkin(9); 0, 0, -1, -pkin(11); t188, t194, 0, 0; 0, 0, 0, 1; t166, -t200, 0, -pkin(14); 0, 0, -1, 0; t200, t166, 0, -pkin(16); 0, 0, 0, 1; t163, -t164, 0, pkin(1); t164, t163, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t167, t169, 0, 0; -t169, t167, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t167, -t169, 0, pkin(2); t169, t167, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; -t176, -t175, 0, cos(pkin(21)) * pkin(4); t175, -t176, 0, -sin(pkin(21)) * pkin(4); 0, 0, 1, 0; 0, 0, 0, 1; t163, t164, 0, t185 * pkin(3); -t164, t163, 0, t182 * pkin(3); 0, 0, 1, 0; 0, 0, 0, 1; t189, -t195, 0, t187 * pkin(6); t195, t189, 0, t184 * pkin(6); 0, 0, 1, 0; 0, 0, 0, 1; t157, t155, 0, t186 * pkin(10); -t155, t157, 0, t183 * pkin(10); 0, 0, 1, 0; 0, 0, 0, 1; 1, 0, 0, pkin(7); 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; 1, 0, 0, pkin(12); 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; -1, 0, 0, pkin(8); 0, -1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,16);             % numerisch
else,                         T_mdh = sym('xx', [4,4,16]); end % symbolisch

for i = 1:16
  T_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
